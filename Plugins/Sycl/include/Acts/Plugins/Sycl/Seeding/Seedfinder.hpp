// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s).
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include <CL/sycl.hpp>

namespace Acts::Sycl {

struct LinCircle {
  float Zo;
  float cotTheta;
  float iDeltaR;
  float Er;
  float U;
  float V;
};

// store SpacePoint data in float arrays, index them with enum values
enum eSpacePoint {
  eX, eY, eZ, eRadius, eVarianceR, eVarianceZ, eSP
};

// store SeedfinderConfig data in float array, index it with enum values
enum eConfigData {
  eDeltaRMin,
  eDeltaRMax,
  eCotThetaMax,
  eCollisionRegionMin,
  eCollisionRegionMax
};

// maximize number of bottom and top spacepoints per middle sp (to be able to load data to global buffer)
enum eMaxData {
  eMaxBottomPerMiddleSP, eMaxTopPerMiddleSP
};

// store linear circle data in float arrays
enum eLinCircle {
  eZo, eCotTheta, eIDeltaR, eEr, eU, eV, eLIN
};

void offloadDupletSearchBottom( const std::vector<float>& configData,
                                const std::vector<int>& maxData,
                                std::vector<int>& indBPerMSpCompat,
                                std::vector<int>& indTPerMSpCompat,
                                std::vector<int>& numBotCompatPerMSP,
                                std::vector<int>& numTopCompatPerMSP,
                                const std::vector<float>& bottomSPs,
                                const std::vector<float>& middleSPs,
                                const std::vector<float>& topSPs);

void offloadTransformCoordiantes( const std::vector<int>& indBPerMSpCompat,
                                  const std::vector<int>& indTPerMSpCompat,
                                  const std::vector<int>& numBotCompatPerMSP,
                                  const std::vector<int>& numTopCompatPerMSP,
                                  const std::vector<float>& bottomSPs,
                                  const std::vector<float>& middleSPs,
                                  const std::vector<float>& topSPs,
                                  std::vector<float>& linCircleBot,
                                  std::vector<float>& linCircleTop);

void outputPlatforms();
void testDevice();

template <typename external_spacepoint_t>
class Seedfinder {
  public:
  Seedfinder(Acts::SeedfinderConfig<external_spacepoint_t> config);

  ~Seedfinder() = default;
  Seedfinder() = delete;
  Seedfinder(const Seedfinder<external_spacepoint_t>&) = delete;
  Seedfinder<external_spacepoint_t>& operator=(
    const Seedfinder<external_spacepoint_t>&) = delete;

  template <typename sp_range_t>
  std::vector<Seed<external_spacepoint_t> > createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const;

 private:

  void transformCoordinates(
    std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    std::vector<LinCircle>& linCircleVec) const;

  Acts::SeedfinderConfig<external_spacepoint_t> m_config;
};

} // namespace Acts

// Include the template implementation.
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.ipp"