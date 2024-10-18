// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
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

/// Subtract bound parameters and take care of angle periodicity for phi and
/// theta.
///
/// @param lhs The left hand side bound parameters
/// @param rhs The right hand side bound parameters
///
/// @return The difference of the bound parameters
inline BoundVector subtractBoundParameters(const BoundVector& lhs,
                                           const BoundVector& rhs) {
  BoundVector result = lhs - rhs;
  result[eBoundPhi] =
      detail::difference_periodic(lhs[eBoundPhi], rhs[eBoundPhi], 2 * M_PI);
  return result;
}

}  // namespace Acts
