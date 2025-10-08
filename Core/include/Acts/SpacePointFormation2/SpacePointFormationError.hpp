// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <system_error>
#include <type_traits>

namespace Acts {

/// Error codes for space point formation operations
enum class SpacePointFormationError {
  // ensure all values are non-zero
  ClusterPairDistanceExceeded = 1,
  ClusterPairThetaDistanceExceeded = 2,
  ClusterPairPhiDistanceExceeded = 3,
  CosmicToleranceNotMet = 4,
  OutsideLimits = 5,
  OutsideRelaxedLimits = 6,
  NoSolutionFound = 7
};

/// Create error code from SpacePointFormationError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(SpacePointFormationError e);

}  // namespace Acts

// register with STL
namespace std {
template <>
struct is_error_code_enum<Acts::SpacePointFormationError> : std::true_type {};
}  // namespace std
