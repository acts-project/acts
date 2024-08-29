// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/SympyStepper.hpp"

#include "Acts/Propagator/detail/SympyCovarianceEngine.hpp"
#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"
#include "Acts/Utilities/QuickMath.hpp"

#include <cmath>
#include <cstdint>

#include "codegen/sympy_stepper_math.hpp"

namespace Acts {

SympyStepper::SympyStepper(std::shared_ptr<const MagneticFieldProvider> bField)
    : m_bField(std::move(bField)) {}

SympyStepper::SympyStepper(const Config& config) : m_bField(config.bField) {}

SympyStepper::State SympyStepper::makeState(
    std::reference_wrapper<const GeometryContext> gctx,
    std::reference_wrapper<const MagneticFieldContext> mctx,
    const BoundTrackParameters& par, double ssize) const {
  return State{gctx, m_bField->makeCache(mctx), par, ssize};
}

void SympyStepper::resetState(State& state, const BoundVector& boundParams,
                              const BoundSquareMatrix& cov,
                              const Surface& surface,
                              const double stepSize) const {
  FreeVector freeParams =
      transformBoundToFreeParameters(surface, state.geoContext, boundParams);

  // Update the stepping state
  state.pars = freeParams;
  state.cov = cov;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

Result<std::tuple<BoundTrackParameters, BoundMatrix, double>>
SympyStepper::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  return detail::sympy::boundState(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated,
      freeToBoundCorrection);
}

std::tuple<CurvilinearTrackParameters, BoundMatrix, double>
SympyStepper::curvilinearState(State& state, bool transportCov) const {
  return detail::sympy::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated);
}

void SympyStepper::update(State& state, const FreeVector& freeParams,
                          const BoundVector& /*boundParams*/,
                          const Covariance& covariance,
                          const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
}

void SympyStepper::update(State& state, const Vector3& uposition,
                          const Vector3& udirection, double qOverP,
                          double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qOverP;
}

void SympyStepper::transportCovarianceToCurvilinear(State& state) const {
  detail::sympy::transportCovarianceToCurvilinear(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars.template segment<3>(eFreeDir0));
}

void SympyStepper::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::sympy::transportCovarianceToBound(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, freeToBoundCorrection);
}

Result<double> SympyStepper::stepImpl(
    State& state, Direction stepDirection, double stepTolerance,
    double stepSizeCutOff, std::size_t maxRungeKuttaStepTrials) const {
  auto pos = position(state);
  auto dir = direction(state);
  double t = time(state);
  double qop = qOverP(state);
  double m = particleHypothesis(state).mass();
  double p_abs = absoluteMomentum(state);

  auto getB = [&](const double* p) -> Result<Vector3> {
    return getField(state, {p[0], p[1], p[2]});
  };

  const auto calcStepSizeScaling = [&](const double errorEstimate_) -> double {
    // For details about these values see ATL-SOFT-PUB-2009-001
    constexpr double lower = 0.25;
    constexpr double upper = 4.0;
    // This is given by the order of the Runge-Kutta method
    constexpr double exponent = 0.25;

    // Whether to use fast power function if available
    constexpr bool tryUseFastPow{false};

    double x = stepTolerance / errorEstimate_;

    if constexpr (exponent == 0.25 && !tryUseFastPow) {
      // This is 3x faster than std::pow
      x = std::sqrt(std::sqrt(x));
    } else if constexpr (std::numeric_limits<double>::is_iec559 &&
                         tryUseFastPow) {
      x = fastPow(x, exponent);
    } else {
      x = std::pow(x, exponent);
    }

    return std::clamp(x, lower, upper);
  };

  double h = state.stepSize.value() * stepDirection;
  double initialH = h;
  std::size_t nStepTrials = 0;
  double errorEstimate = 0.;

  while (true) {
    nStepTrials++;

    // For details about the factor 4 see ATL-SOFT-PUB-2009-001
    Result<bool> res =
        rk4(pos.data(), dir.data(), t, h, qop, m, p_abs, getB, &errorEstimate,
            4 * stepTolerance, state.pars.template segment<3>(eFreePos0).data(),
            state.pars.template segment<3>(eFreeDir0).data(),
            state.pars.template segment<1>(eFreeTime).data(),
            state.derivative.data(),
            state.covTransport ? state.jacTransport.data() : nullptr);
    if (!res.ok()) {
      return res.error();
    }
    // Protect against division by zero
    errorEstimate = std::max(1e-20, errorEstimate);

    if (*res) {
      break;
    }

    const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
    h *= stepSizeScaling;

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(h) < std::abs(stepSizeCutOff)) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
  }

  state.pathAccumulated += h;
  ++state.nSteps;
  state.nStepTrials += nStepTrials;

  const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
  const double nextAccuracy = std::abs(h * stepSizeScaling);
  const double previousAccuracy = std::abs(state.stepSize.accuracy());
  const double initialStepLength = std::abs(initialH);
  if (nextAccuracy < initialStepLength || nextAccuracy > previousAccuracy) {
    state.stepSize.setAccuracy(nextAccuracy);
  }

  return h;
}

void SympyStepper::setIdentityJacobian(State& state) const {
  state.jacobian = BoundMatrix::Identity();
}

}  // namespace Acts
