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

/// Error codes for the multi-space-point track parameter estimation
/// @ingroup errors
enum class TrackParamsEstimationError {
  // ensure all values are non-zero
  /// Fewer than three space points were provided
  NotEnoughSpacePoints = 1,
  /// The fit is degenerate (e.g. all space points coincide)
  DegenerateFit,
};

/// Create error code from @ref TrackParamsEstimationError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(Acts::TrackParamsEstimationError e);

}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::TrackParamsEstimationError> : std::true_type {};
}  // namespace std
