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

/// Recoverable errors in bound-coordinate interpretation.
/// @ingroup errors
enum class BoundParametrizationError {
  /// Inverse direction derivative is singular or insufficiently conditioned.
  SingularDirection = 1,
  /// Canonicalization has no unique local derivative at a pole.
  UndefinedCanonicalization,
};

/// Convert a bound-parametrization error to a standard error code.
/// @param error Error enum value
/// @return Standard error code
std::error_code make_error_code(BoundParametrizationError error);

}  // namespace Acts

namespace std {
template <>
struct is_error_code_enum<Acts::BoundParametrizationError> : std::true_type {};
}  // namespace std
