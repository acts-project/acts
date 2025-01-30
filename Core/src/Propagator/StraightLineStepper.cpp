// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/StraightLineStepper.hpp"

#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"

namespace Acts {

StraightLineStepper::State StraightLineStepper::makeState(
    const Options& options) const {
  State state{options};

  state.stepSize = ConstrainedStep(options.maxStepSize);

  return state;
}

void StraightLineStepper::initialize(State& state,
                                     const BoundTrackParameters& par) const {
  Vector3 position = par.position(state.options.geoContext);
  Vector3 direction = par.direction();

  state.particleHypothesis = par.particleHypothesis();

  state.pars.template segment<3>(eFreePos0) = position;
  state.pars.template segment<3>(eFreeDir0) = direction;
  state.pars[eFreeTime] = par.time();
  state.pars[eFreeQOverP] = par.parameters()[eBoundQOverP];

  // Init the jacobian matrix if needed
  if (par.covariance()) {
    // Get the reference surface for navigation
    const auto& surface = par.referenceSurface();
    // set the covariance transport flag to true and copy
    state.covTransport = true;
    state.cov = BoundSquareMatrix(*par.covariance());
    state.jacToGlobal = surface.boundToFreeJacobian(state.options.geoContext,
                                                    position, direction);
    state.jacobian = BoundMatrix::Identity();
    state.jacTransport = FreeMatrix::Identity();
    state.derivative = FreeVector::Zero();
  }

  state.stepSize = ConstrainedStep(state.options.maxStepSize);

  state.pathAccumulated = 0.;
}

void StraightLineStepper::initialize(State& state,
                                     const BoundVector& boundParams,
                                     const BoundMatrix& cov,
                                     ParticleHypothesis particleHypothesis,
                                     const Surface& surface) const {
  FreeVector freeParams = transformBoundToFreeParameters(
      surface, state.options.geoContext, boundParams);

  state.particleHypothesis = particleHypothesis;

  state.pars = freeParams;

  state.covTransport = true;
  state.cov = cov;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.options.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();

  state.stepSize = ConstrainedStep(state.options.maxStepSize);

  state.pathAccumulated = 0.;
}

Result<std::tuple<BoundTrackParameters, BoundMatrix, double>>
StraightLineStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  return detail::boundState(
      state.options.geoContext, surface, state.cov, state.jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal, state.pars,
      state.particleHypothesis, state.covTransport && transportCov,
      state.pathAccumulated, freeToBoundCorrection);
}

std::tuple<CurvilinearTrackParameters, BoundMatrix, double>
StraightLineStepper::curvilinearState(State& state, bool transportCov) const {
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated);
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
}

void StraightLineStepper::update(State& state, const Vector3& uposition,
                                 const Vector3& udirection, double qop,
                                 double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qop;
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
      state.options.geoContext, surface, state.cov, state.jacobian,
      state.jacTransport, state.derivative, state.jacToGlobal, state.pars,
      freeToBoundCorrection);
}

}  // namespace Acts
