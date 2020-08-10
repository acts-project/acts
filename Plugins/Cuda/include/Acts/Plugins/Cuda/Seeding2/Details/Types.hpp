// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// System include(s).
#include <cstddef>

namespace Acts {
namespace Cuda {
namespace Details {

/// Helper struct describing a spacepoint on the device
struct SpacePoint {
  float x = 0.0f;       ///< x-coordinate in beam system coordinates
  float y = 0.0f;       ///< y-coordinate in beam system coordinates
  float z = 0.0f;       ///< z-coordinate in beam system coordinates
  float radius = 0.0f;  ///< radius in beam system coordinates
  float varianceR = 0.0f;
  float varianceZ = 0.0f;
};  // struct SpacePoint

/// Helper struct summarising the results of the dublet search
struct DubletCounts {
  /// The total number of dublets (M-B and M-T) found
  std::size_t nDublets = 0;
  /// The total number of triplet candidates found
  std::size_t nTriplets = 0;
  /// The maximal number of middle-bottom dublets
  unsigned int maxMBDublets = 0;
  /// The maximal number of middle-top dublets
  unsigned int maxMTDublets = 0;
  /// The maximal number of triplets for any middle SP
  unsigned int maxTriplets = 0;
};  // struct DubletCounts

/// Helper struct holding the linearly transformed coordinates of spacepoints
struct LinCircle {
  float Zo = 0.0f;
  float cotTheta = 0.0f;
  float iDeltaR = 0.0f;
  float Er = 0.0f;
  float U = 0.0f;
  float V = 0.0f;
};  // struct LinCircle

/// Structure used in the CUDA-based triplet finding
struct Triplet {
  std::size_t bottomIndex = static_cast<std::size_t>(-1);
  std::size_t topIndex = static_cast<std::size_t>(-1);
  float impactParameter = 0.0f;
  float invHelixDiameter = 0.0f;
  float weight = 0.0f;
};  // struct Triplet

}  // namespace Details
}  // namespace Cuda
}  // namespace Acts
