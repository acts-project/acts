// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <numeric>
#include <type_traits>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <future>

namespace Acts {

  template <typename external_spacepoint_t, typename platform_t>
  Seedfinder<external_spacepoint_t, platform_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config)
    : m_config(std::move(config)) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;

  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);

  t_metric = std::make_tuple(0,0,0,0);
  
  }
  

  template <typename external_spacepoint_t, typename platform_t>
  template <typename sp_range_t>
  std::vector<Seed<external_spacepoint_t> >  
  Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const{
  std::vector<Seed<external_spacepoint_t>> outputVec;
  auto start_wall = std::chrono::system_clock::now();

  for (auto spM : middleSPs) {    

    auto start_DS = std::chrono::system_clock::now();
    
    float rM = spM->radius();
    float zM = spM->z();
    float varianceRM = spM->varianceR();
    float varianceZM = spM->varianceZ();
    
    // Doublet search    
    auto compatBottomSP =
      SeedfinderCpuFunctions<external_spacepoint_t,
			     sp_range_t>::searchDoublet(true, bottomSPs, *spM, m_config);
    
    // no bottom SP found -> try next spM
    if (compatBottomSP.empty()) {
      continue;
    }

    auto compatTopSP =
      SeedfinderCpuFunctions<external_spacepoint_t,
			     sp_range_t>::searchDoublet(false, topSPs, *spM, m_config);

    // no top SP found -> try next spM
    if (compatTopSP.empty()) {
      continue;
    }
    
    auto end_DS = std::chrono::system_clock::now();
    std::chrono::duration<double> elapse_DS = end_DS-start_DS;
    std::get<0>(t_metric) += elapse_DS.count();
    
    // contains parameters required to calculate circle with linear equation

    auto start_TC = std::chrono::system_clock::now();
    
    // ...for bottom-middle
    std::vector<LinCircle> linCircleBottom;
    // ...for middle-top
    std::vector<LinCircle> linCircleTop;
    
    SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatBottomSP, *spM, true, linCircleBottom);
    SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatTopSP, *spM, false, linCircleTop);

    auto end_TC = std::chrono::system_clock::now();
    std::chrono::duration<double> elapse_TC = end_TC-start_TC;
    std::get<1>(t_metric) += elapse_TC.count();

    auto start_TS = std::chrono::system_clock::now();
    
    auto seedsPerSpM = SeedfinderCpuFunctions<external_spacepoint_t,sp_range_t>::searchTriplet(*spM, compatBottomSP, compatTopSP, linCircleBottom, linCircleTop, m_config);
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);

    auto end_TS = std::chrono::system_clock::now();
    std::chrono::duration<double> elapse_TS = end_TS-start_TS;
    std::get<2>(t_metric) += (elapse_TS).count();
  }  

  auto end_wall = std::chrono::system_clock::now();
  std::chrono::duration<double> elapse_wall = end_wall-start_wall;
  std::get<3>(t_metric) += elapse_wall.count();
  
  return outputVec;
  }

  template< typename external_spacepoint_t, typename platform_t >
  std::tuple< double, double, double, double >
  Seedfinder<external_spacepoint_t, platform_t>::getTimeMetric() { return t_metric; }
  
}// namespace Acts
