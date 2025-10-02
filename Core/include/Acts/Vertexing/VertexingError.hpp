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

/// Error codes for vertexing operations
enum class VertexingError {
  // ensure all values are non-zero
  NumericFailure = 1,
  EmptyInput,
  SeedingError,
  NotConverged,
  ElementNotFound,
  NoCovariance,
  SingularMatrix,
  NonPositiveVariance,
  MatrixNotPositiveDefinite,
  InvalidInput,
  CouldNotRemoveTrack,
};

/// Create error code from VertexingError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(Acts::VertexingError e);

}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::VertexingError> : std::true_type {};
}  // namespace std
