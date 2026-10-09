// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <limits>

namespace Acts {

/// Decode the direction of canonical global-spherical bound coordinates.
/// Accepts vectors, Eigen maps and expressions without copying the full state.
/// @tparam Derived Eigen expression type
/// @param parameters Raw bound coordinates
/// @return Global unit direction vector
/// @note This helper assumes the canonical global-spherical parametrization.
template <typename Derived>
Vector3 directionFromBoundParameters(
    const Eigen::MatrixBase<Derived>& parameters) {
  return makeDirectionFromPhiTheta(parameters[eBoundDir0],
                                   parameters[eBoundDir1]);
}

/// Check if a bound vector is valid. This checks the following:
/// - All values are finite
/// - (optionally) The phi value is in the range [-pi, pi)
/// - (optionally) The theta value is in the range [0, pi]
///
/// @param v The bound vector to check
/// @param validateAngleRange If true, the phi and theta values are range checked
/// @param epsilon The epsilon to use for the checks
/// @param maxAbsEta The maximum allowed eta value
///
/// @return True if the bound vector is valid
bool isBoundVectorValid(
    const BoundVector& v, bool validateAngleRange, double epsilon = 1e-6,
    double maxAbsEta = std::numeric_limits<double>::infinity());

/// Check if a free vector is valid. This checks the following:
/// - All values are finite
/// - Direction is normalized
///
/// @param v The free vector to check
/// @param epsilon The epsilon to use for the checks
/// @param maxAbsEta The maximum allowed eta value
///
/// @return True if the free vector is valid
bool isFreeVectorValid(
    const FreeVector& v, double epsilon = 1e-6,
    double maxAbsEta = std::numeric_limits<double>::infinity());

/// Normalize the canonical global-spherical bound parameter angles
///
/// @param boundParams The bound parameters to normalize
///
/// @return The normalized bound parameters
inline BoundVector normalizeBoundParameters(const BoundVector& boundParams) {
  BoundVector result = boundParams;
  std::tie(result[eBoundDir0], result[eBoundDir1]) =
      detail::normalizePhiTheta(result[eBoundDir0], result[eBoundDir1]);
  return result;
}

/// Add canonical global-spherical bound parameters, respecting angle
/// periodicity. This is intended for small differences only i.e. KF updates.
///
/// @param lhs The left hand side bound parameters
/// @param rhs The right hand side bound parameters
///
/// @return The sum of the bound parameters
inline BoundVector addBoundParameters(const BoundVector& lhs,
                                      const BoundVector& rhs) {
  return normalizeBoundParameters(lhs + rhs);
}

/// Subtract canonical global-spherical bound parameters, respecting angle
/// periodicity. This is intended for small differences only i.e. KF updates.
///
/// @param lhs The left hand side bound parameters
/// @param rhs The right hand side bound parameters
///
/// @return The difference of the bound parameters
inline BoundVector subtractBoundParameters(const BoundVector& lhs,
                                           const BoundVector& rhs) {
  BoundVector result = lhs - rhs;
  result[eBoundDir0] = detail::radian_sym(result[eBoundDir0]);
  result[eBoundDir1] = detail::radian_sym(result[eBoundDir1]);
  return result;
}

}  // namespace Acts
