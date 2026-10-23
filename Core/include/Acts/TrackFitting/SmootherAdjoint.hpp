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

namespace Acts {

/// @addtogroup track_fitting
/// @{

/// The information a smoother carries out of its backward walk.
///
/// After the walk the adjoint is transported through the jacobian of the last
/// smoothed state, so it lives on the node that state was predicted from. For
/// a complete track that node is the start parameters; for a track segment
/// that was filtered from the state of another track it is that shared state.
///
/// The notation follows the modified Bryson-Frazier smoother: a node with the
/// estimate `(x, P)` and no further measurements of its own is smoothed as
/// `x' = x - P smallLambda` and `P' = P - P bigLambda P`.
struct SmootherAdjoint {
  /// Information matrix term of the adjoint
  BoundMatrix bigLambda = BoundMatrix::Zero();
  /// Information vector term of the adjoint
  BoundVector smallLambda = BoundVector::Zero();
};

/// Apply a carried adjoint to the estimate of the node it was transported to.
///
/// @param mean The estimate at the node before the smoothed segment
/// @param cov Its covariance
/// @param adjoint The adjoint carried out of the segment
/// @return The estimate at the node with the segment's measurements included
inline std::pair<BoundVector, BoundMatrix> applySmootherAdjoint(
    const BoundVector& mean, const BoundMatrix& cov,
    const SmootherAdjoint& adjoint) {
  const BoundVector newMean =
      normalizeBoundParameters(mean - cov * adjoint.smallLambda);
  const BoundMatrix newCov = cov - cov * adjoint.bigLambda * cov;
  // keep the covariance symmetric against rounding
  return {newMean, 0.5 * (newCov + newCov.transpose())};
}

/// @}

}  // namespace Acts
