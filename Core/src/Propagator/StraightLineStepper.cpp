// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/StraightLineStepper.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"

namespace Acts {

StraightLineStepper::State StraightLineStepper::makeState(
    const Options& options) const {
  State state{options};
  return state;
}

void StraightLineStepper::initialize(State& state,
                                     const BoundParameters& par) const {
  initialize(state, par.parameters(), par.covariance(),
             par.particleHypothesis(), par.referenceSurface());
}

void StraightLineStepper::initialize(State& state,
                                     const BoundVector& boundParams,
                                     const std::optional<BoundMatrix>& cov,
                                     ParticleHypothesis particleHypothesis,
                                     const Surface& surface) const {
  FreeVector freeParams = transformBoundToFreeParameters(
      surface, state.options.geoContext, boundParams);

  state.particleHypothesis = particleHypothesis;

  state.pathAccumulated = 0;
  state.nSteps = 0;
  state.stepSize = ConstrainedStep();
  state.stepSize.setAccuracy(state.options.initialStepSize);
  state.stepSize.setUser(state.options.maxStepSize);
  state.previousStepSize = 0;
  state.statistics = StepperStatistics();

  state.pars = freeParams;

  // Init the jacobian matrix if needed
  state.covTransport = cov.has_value();
  if (state.covTransport) {
    state.cov = *cov;
    state.jacToGlobal = surface.boundToFreeJacobian(
        state.options.geoContext, freeParams.segment<3>(eFreePos0),
        freeParams.segment<3>(eFreeDir0));
    state.jacTransport = FreeMatrix::Identity();
    state.derivative = FreeVector::Zero();
  }
}

Result<StraightLineStepper::BoundParameters>
StraightLineStepper::boundParameters(const State& state,
                                     const Surface& surface) const {
  return detail::boundParameters(
      state.options.geoContext, surface, state.pars,
      state.covTransport ? std::optional(state.cov) : std::nullopt,
      state.particleHypothesis);
}

StraightLineStepper::BoundParameters StraightLineStepper::curvilinearParameters(
    const State& state) const {
  return detail::curvilinearParameters(
      state.pars, state.covTransport ? std::optional(state.cov) : std::nullopt,
      state.particleHypothesis);
}

void StraightLineStepper::update(State& state, const FreeVector& freeParams,
                                 const BoundVector& /*boundParams*/,
                                 const Covariance& covariance,
                                 const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.options.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

void StraightLineStepper::update(State& state, const Vector3& uposition,
                                 const Vector3& udirection, double qop,
                                 double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qop;
}

StraightLineStepper::Jacobian StraightLineStepper::transportToCurvilinear(
    State& state) const {
  Jacobian jacobian = Jacobian::Identity();
  if (!state.covTransport) {
    return jacobian;
  }
  detail::transportCovarianceToCurvilinear(
      state.cov, jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, std::nullopt,
      state.pars.template segment<3>(eFreeDir0));
  return jacobian;
}

Result<StraightLineStepper::Jacobian> StraightLineStepper::transportToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  if (!surface.isOnSurface(state.options.geoContext, position(state),
                           direction(state), BoundaryTolerance::Infinite())) {
    return Result<Jacobian>::failure(SurfaceError::GlobalPositionNotOnSurface);
  }

  Jacobian jacobian = Jacobian::Identity();
  if (!state.covTransport) {
    return Result<Jacobian>::success(jacobian);
  }
  detail::transportCovarianceToBound(
      state.options.geoContext, surface, state.cov, jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal, std::nullopt,
      state.pars, freeToBoundCorrection);
  return Result<Jacobian>::success(jacobian);
}

}  // namespace Acts
