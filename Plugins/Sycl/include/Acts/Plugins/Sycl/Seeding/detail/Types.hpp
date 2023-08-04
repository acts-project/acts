// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// System include(s)
#include <cstddef>
#include <cstdint>

namespace Acts::Sycl::detail {
/// Simplified internal space point object with parameters required for the
/// seeding algorithm on the device side.
struct DeviceSpacePoint {
  float x;
  float y;
  float z;
  /// radius
  float r;
  /// variance R
  float varR;
  /// variance Z
  float varZ;
};

/// Simple structure holding information about duplet circles (middle-bottom and
/// middle-top) in linear equation form.
/// They are used on the device side in SYCL kernels.
struct DeviceLinEqCircle {
  float zo;
  float cotTheta;
  float iDeltaR;
  float er;
  float u;
  float v;
};

/// Config parameters of SeedFinder and SeedFilter classes
/// needed for the seeding algorithm on the device side.
struct DeviceSeedFinderConfig {
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

/// Struct holding information about triplets that are calculated in the triplet
/// search kernel and needed in the next kernel (triplet filter).
struct DeviceTriplet {
  float curvature;
  float impact;
  int32_t topSPIndex;
};

/// Struct returned from the SYCL seeding algorithm to the host
///
/// Stores the indices of the space points and the weight of the seed.
/// They index arrays on the host side, constructed in
/// @c Acts::Sycl::SeedFinder::createSeedsForGroup(...).
struct SeedData {
  uint32_t bottom;
  uint32_t top;
  uint32_t middle;
  float weight;
};
}  // namespace Acts::Sycl::detail
