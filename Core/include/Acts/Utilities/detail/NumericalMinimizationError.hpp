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

namespace Acts::detail {

/// Error codes for the numerical minimization routines in
/// @c NumericalMinimization.hpp
/// @ingroup errors
enum class NumericalMinimizationError {
  // ensure all values are non-zero
  /// The Nelder-Mead simplex search did not converge within the iteration cap
  DidNotConverge = 1,
  /// No vertex of the initial simplex evaluated to a finite objective value
  NoFiniteVertex,
  /// A finite-difference step was not strictly positive
  InvalidStep,
  /// The finite-difference Hessian was not positive definite, i.e. the probed
  /// point is not a genuine minimum
  NotPositiveDefinite,
};

/// Create error code from NumericalMinimizationError
/// @param e The error code enum value
/// @return Standard error code
std::error_code make_error_code(NumericalMinimizationError e);

}  // namespace Acts::detail

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::detail::NumericalMinimizationError>
    : std::true_type {};
}  // namespace std
