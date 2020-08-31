// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"

#include <algorithm>
#include <cmath>

namespace Acts::Sycl {
using namespace Acts::Sycl::detail;
template <typename external_spacepoint_t>
Seedfinder<external_spacepoint_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config,
    const Acts::Sycl::DeviceExperimentCuts& cuts,
    const std::string& device_name_substring)
    : m_config(config), m_deviceCuts(cuts) {
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

  m_queue = createQueue(device_name_substring);
}

template <typename external_spacepoint_t>
template <typename sp_range_t>
std::vector<Acts::Seed<external_spacepoint_t>>
Seedfinder<external_spacepoint_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  std::vector<deviceSpacePoint> deviceBottomSPs;
  std::vector<deviceSpacePoint> deviceMiddleSPs;
  std::vector<deviceSpacePoint> deviceTopSPs;

  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>
      bottomSPvec;
  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>
      middleSPvec;
  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*> topSPvec;

  for (auto SP : bottomSPs) {
    bottomSPvec.push_back(SP);
    deviceBottomSPs.insert(
        deviceBottomSPs.end(),
        deviceSpacePoint{SP->x(), SP->y(), SP->z(), SP->radius(),
                         SP->varianceR(), SP->varianceZ()});
  }

  for (auto SP : middleSPs) {
    middleSPvec.push_back(SP);
    deviceMiddleSPs.insert(
        deviceMiddleSPs.end(),
        deviceSpacePoint{SP->x(), SP->y(), SP->z(), SP->radius(),
                         SP->varianceR(), SP->varianceZ()});
  }

  for (auto SP : topSPs) {
    topSPvec.push_back(SP);
    deviceTopSPs.insert(
        deviceTopSPs.end(),
        deviceSpacePoint{SP->x(), SP->y(), SP->z(), SP->radius(),
                         SP->varianceR(), SP->varianceZ()});
  }

  const deviceSeedfinderConfig deviceConfigData = {
      m_config.deltaRMin,
      m_config.deltaRMax,
      m_config.cotThetaMax,
      m_config.collisionRegionMin,
      m_config.collisionRegionMax,
      m_config.maxScatteringAngle2,
      m_config.sigmaScattering,
      m_config.minHelixDiameter2,
      m_config.pT2perRadius,
      m_config.seedFilter->getSeedFilterConfig().deltaInvHelixDiameter,
      m_config.seedFilter->getSeedFilterConfig().impactWeightFactor,
      m_config.seedFilter->getSeedFilterConfig().deltaRMin,
      m_config.seedFilter->getSeedFilterConfig().compatSeedWeight,
      m_config.impactMax,
      m_config.seedFilter->getSeedFilterConfig().compatSeedLimit,
  };

  std::vector<std::vector<SeedData>> seeds;

  createSeedsForGroupSycl(m_queue, deviceConfigData, m_deviceCuts,
                          deviceBottomSPs, deviceMiddleSPs, deviceTopSPs,
                          seeds);

  for (size_t mi = 0; mi < seeds.size(); ++mi) {
    std::vector<std::pair<
        float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
        seedsPerSPM;
    for (size_t j = 0; j < seeds[mi].size(); ++j) {
      auto& bottomSP = *(bottomSPvec[seeds[mi][j].bottom]);
      auto& middleSP = *(middleSPvec[mi]);
      auto& topSP = *(topSPvec[seeds[mi][j].top]);
      float weight = seeds[mi][j].weight;

      seedsPerSPM.emplace_back(std::make_pair(
          weight, std::make_unique<const InternalSeed<external_spacepoint_t>>(
                      bottomSP, middleSP, topSP, 0)));
    }
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSPM, outputVec);
  }
  return outputVec;
}
}  // namespace Acts::Sycl