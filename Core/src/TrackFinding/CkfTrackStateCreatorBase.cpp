// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/CkfTrackStateCreatorBase.hpp"

#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <cmath>

namespace Acts::detail {

Result<CkfTypes::BoundState> transportBoundStateStraightLine(
    const GeometryContext& gctx, const CkfTypes::BoundState& boundState,
    const Surface& surface, Direction propagationDirection) {
  const auto& [params, jacobian, pathLength] = boundState;
  if (!params.covariance().has_value()) {
    return PropagatorError::Failure;
  }

  const Surface& startSurface = params.referenceSurface();
  const FreeVector startParams =
      transformBoundToFreeParameters(startSurface, gctx, params.parameters());
  const Vector3 position = startParams.segment<3>(eFreePos0);
  const Vector3 direction = startParams.segment<3>(eFreeDir0);

  // the target has to be ahead in the propagation direction
  const auto intersection =
      surface
          .intersect(gctx, position, propagationDirection * direction,
                     BoundaryTolerance::Infinite())
          .closestForward();
  if (!intersection.isValid()) {
    return PropagatorError::Failure;
  }
  // the path length along the momentum, as the stepper counts it
  const double s = propagationDirection * intersection.pathLength();

  // dt/ds and its derivative with respect to q/p
  const double qOverP = params.qOverP();
  const double mass = params.particleHypothesis().mass();
  double dtds = 1;
  double dtdsDQOverP = 0;
  if (qOverP != 0) {
    const double p = params.particleHypothesis().extractMomentum(qOverP);
    const double massOverP2 = (mass * mass) / (p * p);
    dtds = std::sqrt(1 + massOverP2);
    dtdsDQOverP = massOverP2 / (qOverP * dtds);
  }

  FreeVector endParams = startParams;
  endParams.segment<3>(eFreePos0) = position + s * direction;
  endParams[eFreeTime] += s * dtds;

  FreeMatrix freeTransportJacobian = FreeMatrix::Identity();
  freeTransportJacobian.block<3, 3>(eFreePos0, eFreeDir0) =
      s * SquareMatrix<3>::Identity();
  freeTransportJacobian(eFreeTime, eFreeQOverP) = s * dtdsDQOverP;

  FreeVector freeToPathDerivatives = FreeVector::Zero();
  freeToPathDerivatives.segment<3>(eFreePos0) = direction;
  freeToPathDerivatives[eFreeTime] = dtds;

  const BoundToFreeMatrix boundToFreeJacobian =
      startSurface.boundToFreeJacobian(gctx, position, direction);
  const BoundMatrix transportJacobian = boundToBoundTransportJacobian(
      gctx, endParams, boundToFreeJacobian, freeTransportJacobian,
      freeToPathDerivatives, surface);

  Result<BoundVector> endBoundParams =
      transformFreeToBoundParameters(endParams, surface, gctx);
  if (!endBoundParams.ok()) {
    return endBoundParams.error();
  }

  const BoundMatrix covariance = transportJacobian * (*params.covariance()) *
                                 transportJacobian.transpose();

  return CkfTypes::BoundState{
      BoundTrackParameters(surface.getSharedPtr(), *endBoundParams, covariance,
                           params.particleHypothesis()),
      transportJacobian * jacobian, pathLength + s};
}

}  // namespace Acts::detail
