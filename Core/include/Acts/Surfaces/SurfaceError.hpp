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

/// Error codes for surface operations
/// @ingroup errors
enum class SurfaceError {
  // ensure all values are non-zero
  /// Global to local transformation failed: position not on surface.
  GlobalPositionNotOnSurface = 1,
};

/// Create error code from SurfaceError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(Acts::SurfaceError e);

}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::SurfaceError> : std::true_type {};
}  // namespace std
