// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

namespace Acts {

/// Normalize the bound parameter angles
///
/// @param boundParams The bound parameters to normalize
///
/// @return The normalized bound parameters
inline BoundVector normalizeBoundParameters(const BoundVector& boundParams) {
  BoundVector result = boundParams;
  std::tie(result[eBoundPhi], result[eBoundTheta]) =
      detail::normalizePhiTheta(result[eBoundPhi], result[eBoundTheta]);
  return result;
}

/// Add bound parameters and take care of angle periodicity for phi and theta.
/// This is intended for small differences only i.e. KF updates.
///
/// @param lhs The left hand side bound parameters
/// @param rhs The right hand side bound parameters
///
/// @return The sum of the bound parameters
inline BoundVector addBoundParameters(const BoundVector& lhs,
                                      const BoundVector& rhs) {
  return normalizeBoundParameters(lhs + rhs);
}

/// Subtract bound parameters and take care of angle periodicity for phi and
/// theta. This is intended for small differences only i.e. KF updates.
///
/// @param lhs The left hand side bound parameters
/// @param rhs The right hand side bound parameters
///
/// @return The difference of the bound parameters
inline BoundVector subtractBoundParameters(const BoundVector& lhs,
                                           const BoundVector& rhs) {
  BoundVector result = lhs - rhs;
  result[eBoundPhi] = detail::radian_sym(result[eBoundPhi]);
  result[eBoundTheta] = detail::radian_sym(result[eBoundTheta]);
  return result;
}

}  // namespace Acts
