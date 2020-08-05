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

// store SpacePoint data in float arrays, index them with enum values
enum eSpacePoint {
  eX = 0,
  eY = 1,
  eZ = 2,
  eRadius = 3,
  eVarianceR = 4,
  eVarianceZ = 5,
  eSP = 6
};

// store SeedfinderConfig data in float array, index it with enum values
enum eConfigData {
  eDeltaRMin = 0,
  eDeltaRMax = 1,
  eCotThetaMax = 2,
  eCollisionRegionMin = 3,
  eCollisionRegionMax = 4,
  eMaxScatteringAngle2 = 5,
  eSigmaScattering = 6,
  eMinHelixDiameter2 = 7,
  ePT2perRadius = 8,
  eDeltaInvHelixDiameter = 9,
  eImpactWeightFactor = 10,
  eFilterDeltaRMin = 11,
  eCompatSeedWeight = 12,
  eCompatSeedLimit = 13
};

// store linear circle data in float arrays, index them with enum values
enum eLinCircle {
  eZo = 0,
  eCotTheta = 1,
  eIDeltaR = 2,
  eEr = 3,
  eU = 4,
  eV = 5,
  eLIN = 6
};

void offloadComputations( cl::sycl::queue q,
                          const std::vector<float>& configData,
                          const std::vector<float>& bottomSPs,
                          const std::vector<float>& middleSPs,
                          const std::vector<float>& topSPs,
                          std::vector<std::vector<int>>& seedIndices,
                          std::vector<float>& seedWeight);

void outputPlatforms();
void testDevice();

struct nvidia_selector : public cl::sycl::device_selector {
  int operator()(const cl::sycl::device& d) const override {
    if(d.get_info<cl::sycl::info::device::vendor>().find("NVIDIA") != std::string::npos) {
      return 1;
    }
    else {
      return -1;
    }
  }; 
};

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

  Acts::SeedfinderConfig<external_spacepoint_t> m_config;
  cl::sycl::queue m_queue;
};

} // namespace Acts::Sycl

// Include the template implementation.
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.ipp"