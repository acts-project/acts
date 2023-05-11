// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/StraightLineStepper.hpp"

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"

namespace Acts {

Result<std::tuple<BoundTrackParameters, BoundMatrix, double>>
StraightLineStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  return detail::boundState(
      state.geoContext, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars,
      state.covTransport and transportCov, state.pathAccumulated, surface,
      freeToBoundCorrection);
}

std::tuple<CurvilinearTrackParameters, BoundMatrix, double>
StraightLineStepper::curvilinearState(State& state, bool transportCov) const {
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.covTransport and transportCov,
      state.pathAccumulated);
}

void StraightLineStepper::update(State& state, const FreeVector& freeParams,
                                 const BoundVector& boundParams,
                                 const Covariance& covariance,
                                 const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal =
      surface.boundToFreeJacobian(state.geoContext, boundParams);
}

void StraightLineStepper::update(State& state, const Vector3& uposition,
                                 const Vector3& udirection, double up,
                                 double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = (state.q != 0. ? state.q / up : 1. / up);
}

void StraightLineStepper::transportCovarianceToCurvilinear(State& state) const {
  detail::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars.template segment<3>(eFreeDir0));
}

void StraightLineStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::transportCovarianceToBound(
      state.geoContext, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, surface,
      freeToBoundCorrection);
}

void StraightLineStepper::resetState(State& state,
                                     const BoundVector& boundParams,
                                     const BoundSymMatrix& cov,
                                     const Surface& surface,
                                     const Direction navDir,
                                     const double stepSize) const {
  // Update the stepping state
  update(state,
         detail::transformBoundToFreeParameters(surface, state.geoContext,
                                                boundParams),
         boundParams, cov, surface);
  state.navDir = navDir;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal =
      surface.boundToFreeJacobian(state.geoContext, boundParams);
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

}  // namespace Acts
