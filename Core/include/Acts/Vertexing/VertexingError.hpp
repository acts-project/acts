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
/// @ingroup errors
enum class VertexingError {
  // ensure all values are non-zero
  /// Numeric failure in calculation.
  NumericFailure = 1,
  /// Empty input provided.
  EmptyInput,
  /// Error while finding vertex seed.
  SeedingError,
  /// Unable to converge.
  NotConverged,
  /// Unable to find element.
  ElementNotFound,
  /// No covariance provided.
  NoCovariance,
  /// Encountered non-invertible matrix.
  SingularMatrix,
  /// Encountered negative or zero variance.
  NonPositiveVariance,
  /// Encountered a matrix that is not positive definite.
  MatrixNotPositiveDefinite,
  /// Invalid input provided.
  InvalidInput,
  /// Could not remove track from collection.
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
