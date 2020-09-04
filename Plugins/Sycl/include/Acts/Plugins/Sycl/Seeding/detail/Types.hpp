// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cstddef>
#include <cstdint>

namespace Acts::Sycl::detail {
struct DeviceSpacePoint {
  float x;
  float y;
  float z;
  float r;
  float varR;
  float varZ;
};

// Parameters required to calculate circle with linear equation.
struct DeviceLinEqCircle {
  float zo;
  float cotTheta;
  float iDeltaR;
  float er;
  float u;
  float v;
};

// Predefined parameters for Seedfinder and SeedFilter classes
struct DeviceSeedfinderConfig {
  float deltaRMin;
  float deltaRMax;
  float cotThetaMax;
  float collisionRegionMin;
  float collisionRegionMax;
  float maxScatteringAngle2;
  float sigmaScattering;
  float minHelixDiameter2;
  float pT2perRadius;
  float deltaInvHelixDiameter;
  float impactWeightFactor;
  float filterDeltaRMin;
  float compatSeedWeight;
  float impactMax;
  size_t compatSeedLimit;
};

// For SYCL
struct TripletData {
  float curvature;
  float impact;
};

struct SeedData {
  uint32_t bottom;
  uint32_t top;
  uint32_t middle;
  float weight;
};
}  // namespace Acts::Sycl::detail
