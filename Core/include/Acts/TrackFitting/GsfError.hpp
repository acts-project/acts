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

/// Error codes for Gaussian Sum Filter operations
/// @ingroup errors
enum class GsfError {
  /// Start parameters have no Covariance
  StartParametersHaveNoCovariance,
  /// No measurement states found in the forward pass
  NoMeasurementStatesCreatedForward,
  /// No measurement states found in the backward pass
  NoMeasurementStatesCreatedBackward,
  /// No measurement states in the final trajectory
  NoMeasurementStatesCreatedFinal,
};

/// Create error code from @ref GsfError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(GsfError e);

}  // namespace Acts

// register with STL
namespace std {
template <>
struct is_error_code_enum<Acts::GsfError> : std::true_type {};
}  // namespace std
