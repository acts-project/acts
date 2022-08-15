// This file is part of the Acts project.
//
// Copyright (C) 2019-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"

template <typename E, typename A>
Acts::EigenStepper<E, A>::EigenStepper(
    std::shared_ptr<const MagneticFieldProvider> bField)
    : m_bField(std::move(bField)) {}

template <typename E, typename A>
template <typename charge_t>
auto Acts::EigenStepper<E, A>::makeState(
    std::reference_wrapper<const GeometryContext> gctx,
    std::reference_wrapper<const MagneticFieldContext> mctx,
    const SingleBoundTrackParameters<charge_t>& par, NavigationDirection ndir,
    double ssize, double stolerance) const -> State {
  return State{gctx, m_bField->makeCache(mctx), par, ndir, ssize, stolerance};
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::resetState(State& state,
                                          const BoundVector& boundParams,
                                          const BoundSymMatrix& cov,
                                          const Surface& surface,
                                          const NavigationDirection navDir,
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

template <typename E, typename A>
auto Acts::EigenStepper<E, A>::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const
    -> Result<BoundState> {
  return detail::boundState(
      state.geoContext, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars,
      state.covTransport && transportCov, state.pathAccumulated, surface,
      freeToBoundCorrection);
}

template <typename E, typename A>
auto Acts::EigenStepper<E, A>::curvilinearState(State& state,
                                                bool transportCov) const
    -> CurvilinearState {
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.covTransport && transportCov,
      state.pathAccumulated);
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::update(State& state,
                                      const FreeVector& freeParams,
                                      const BoundVector& boundParams,
                                      const Covariance& covariance,
                                      const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal =
      surface.boundToFreeJacobian(state.geoContext, boundParams);
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::update(State& state, const Vector3& uposition,
                                      const Vector3& udirection, double up,
                                      double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = (state.q != 0. ? state.q / up : 1. / up);
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::transportCovarianceToCurvilinear(
    State& state) const {
  detail::transportCovarianceToCurvilinear(state.cov, state.jacobian,
                                           state.jacTransport, state.derivative,
                                           state.jacToGlobal, direction(state));
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::transportCovarianceToBound(
      state.geoContext.get(), state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, surface,
      freeToBoundCorrection);
}

template <typename E, typename A>
template <typename propagator_state_t>
Acts::Result<double> Acts::EigenStepper<E, A>::step(
    propagator_state_t& state) const {
  using namespace UnitLiterals;

  // Runge-Kutta integrator state
  auto& sd = state.stepping.stepData;
  double error_estimate = 0.;
  double h2, half_h;

  auto pos = position(state.stepping);
  auto dir = direction(state.stepping);

  // First Runge-Kutta point (at current position)
  auto fieldRes = getField(state.stepping, pos);
  if (!fieldRes.ok()) {
    return fieldRes.error();
  }
  sd.B_first = *fieldRes;
  if (!state.stepping.extension.validExtensionForStep(state, *this) ||
      !state.stepping.extension.k1(state, *this, sd.k1, sd.B_first, sd.kQoP)) {
    return 0.;
  }

  // The following functor starts to perform a Runge-Kutta step of a certain
  // size, going up to the point where it can return an estimate of the local
  // integration error. The results are stated in the local variables above,
  // allowing integration to continue once the error is deemed satisfactory
  const auto tryRungeKuttaStep =
      [&](const ConstrainedStep& step) -> Result<bool> {
    // helpers because bool and std::error_code are ambiguous
    constexpr auto success = &Result<bool>::success;
    constexpr auto failure = &Result<bool>::failure;

    const double h = step.value();
    // State the square and half of the step size
    h2 = h * h;
    half_h = h * 0.5;

    // Second Runge-Kutta point
    const Vector3 pos1 = pos + half_h * dir + h2 * 0.125 * sd.k1;
    auto field = getField(state.stepping, pos1);
    if (!field.ok()) {
      return failure(field.error());
    }
    sd.B_middle = *field;

    if (!state.stepping.extension.k2(state, *this, sd.k2, sd.B_middle, sd.kQoP,
                                     half_h, sd.k1)) {
      return success(false);
    }

    // Third Runge-Kutta point
    if (!state.stepping.extension.k3(state, *this, sd.k3, sd.B_middle, sd.kQoP,
                                     half_h, sd.k2)) {
      return success(false);
    }

    // Last Runge-Kutta point
    const Vector3 pos2 = pos + h * dir + h2 * 0.5 * sd.k3;
    field = getField(state.stepping, pos2);
    if (!field.ok()) {
      return failure(field.error());
    }
    sd.B_last = *field;
    if (!state.stepping.extension.k4(state, *this, sd.k4, sd.B_last, sd.kQoP, h,
                                     sd.k3)) {
      return success(false);
    }

    // Compute and check the local integration error estimate
    error_estimate =
        h2 * ((sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>() +
              std::abs(sd.kQoP[0] - sd.kQoP[1] - sd.kQoP[2] + sd.kQoP[3]));
    error_estimate = std::max(error_estimate, 1e-20);

    return success(error_estimate <= state.options.tolerance);
  };

  double stepSizeScaling = 1.;
  size_t nStepTrials = 0;
  // Select and adjust the appropriate Runge-Kutta step size as given
  // ATL-SOFT-PUB-2009-001
  while (true) {
    auto res = tryRungeKuttaStep(state.stepping.stepSize);
    if (!res.ok()) {
      return res.error();
    }
    if (!!res.value()) {
      break;
    }

    stepSizeScaling =
        std::min(std::max(0.25f, std::sqrt(std::sqrt(static_cast<float>(
                                     state.options.tolerance /
                                     std::abs(2. * error_estimate))))),
                 4.0f);
    state.stepping.stepSize.scale(stepSizeScaling);

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(state.stepping.stepSize.value()) <
        std::abs(state.options.stepSizeCutOff)) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > state.options.maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
    nStepTrials++;
  }

  // use the adjusted step size
  const double h = state.stepping.stepSize.value();

  // When doing error propagation, update the associated Jacobian matrix
  if (state.stepping.covTransport) {
    // The step transport matrix in global coordinates
    FreeMatrix D;
    if (!state.stepping.extension.finalize(state, *this, h, D)) {
      return EigenStepperError::StepInvalid;
    }

    // for moment, only update the transport part
    state.stepping.jacTransport = D * state.stepping.jacTransport;
  } else {
    if (!state.stepping.extension.finalize(state, *this, h)) {
      return EigenStepperError::StepInvalid;
    }
  }

  // Update the track parameters according to the equations of motion
  state.stepping.pars.template segment<3>(eFreePos0) +=
      h * dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
  state.stepping.pars.template segment<3>(eFreeDir0) +=
      h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
  (state.stepping.pars.template segment<3>(eFreeDir0)).normalize();

  if (state.stepping.covTransport) {
    state.stepping.derivative.template head<3>() =
        state.stepping.pars.template segment<3>(eFreeDir0);
    state.stepping.derivative.template segment<3>(4) = sd.k4;
  }
  state.stepping.pathAccumulated += h;
  if (state.stepping.stepSize.currentType() ==
      ConstrainedStep::Type::accuracy) {
    stepSizeScaling = std::min(
        std::max(0.25f,
                 std::sqrt(std::sqrt(static_cast<float>(
                     state.options.tolerance / std::abs(error_estimate))))),
        4.0f);
    state.stepping.stepSize.scale(stepSizeScaling);
  }

  state.stepping.stepSize.nStepTrials = nStepTrials;

  return h;
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::setIdentityJacobian(State& state) const {
  state.jacobian = BoundMatrix::Identity();
}
