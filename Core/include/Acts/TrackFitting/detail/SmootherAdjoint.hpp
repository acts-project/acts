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
#include "Acts/EventData/TrackParameterHelpers.hpp"

#include <utility>

namespace Acts::detail {

/// The adjoint a smoother carries out of its backward walk, transported onto
/// the node the last smoothed state was predicted from. In the notation of the
/// modified Bryson-Frazier smoother that node is smoothed as
/// `x' = x - P smallLambda`, `P' = P - P bigLambda P`.
struct SmootherAdjoint {
  /// Information matrix term
  BoundMatrix bigLambda = BoundMatrix::Zero();
  /// Information vector term
  BoundVector smallLambda = BoundVector::Zero();
};

/// Apply a carried adjoint to the estimate of the node it was transported to.
///
/// @param mean The estimate at the node
/// @param cov Its covariance
/// @param adjoint The adjoint
/// @return The estimate with the smoothed segment included
inline std::pair<BoundVector, BoundMatrix> applySmootherAdjoint(
    const BoundVector& mean, const BoundMatrix& cov,
    const SmootherAdjoint& adjoint) {
  const BoundVector newMean =
      normalizeBoundParameters(mean - cov * adjoint.smallLambda);
  const BoundMatrix newCov = cov - cov * adjoint.bigLambda * cov;
  // keep the covariance symmetric against rounding
  return {newMean, 0.5 * (newCov + newCov.transpose())};
}

}  // namespace Acts::detail
