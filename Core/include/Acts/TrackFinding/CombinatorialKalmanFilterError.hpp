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

/// Error codes for combinatorial Kalman filter operations
/// @ingroup errors
enum class CombinatorialKalmanFilterError {
  // ensure all values are non-zero
  /// Kalman update failed
  UpdateFailed = 1,
  /// Kalman smooth failed
  SmoothFailed,
  /// Kalman output conversion failed
  OutputConversionFailed,
  /// Measurement selection failed
  MeasurementSelectionFailed,
  /// Propagation reaches max steps before track finding is finished
  PropagationReachesMaxSteps,
  /// No measurement expected on the current surface
  NoMeasurementExpected
};

/// Create error code from @ref CombinatorialKalmanFilterError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(Acts::CombinatorialKalmanFilterError e);

}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::CombinatorialKalmanFilterError>
    : std::true_type {};
}  // namespace std
