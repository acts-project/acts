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

/// Error codes for navigator operations
/// @ingroup errors
enum class NavigatorError {
  // ensure all values are non-zero
  /// We did not end up inside the volume.
  NotInsideExpectedVolume = 1,
  /// Stepper not on surface
  NotOnExpectedSurface = 2,
  /// No start volume could be resolved
  NoStartVolume = 3,
};

/// Create error code from NavigatorError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(Acts::NavigatorError e);

}  // namespace Acts

// register with STL
namespace std {
template <>
struct is_error_code_enum<Acts::NavigatorError> : std::true_type {};
}  // namespace std
