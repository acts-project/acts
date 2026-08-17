// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// The step body, instantiated once per mode.

#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Propagator/detail/SympyStepperDenseStep.hpp"
#include "Acts/Propagator/detail/SympyStepperStatus.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <span>

#include "SympyStepperStep.hpp"
#include "codegen/sympy_stepper_math.hpp"

namespace Acts::detail {

template <SympyStepMode Mode>
Result<double> sympyStep(const SympyStepper& stepper,
                         SympyStepper::State& state, Direction propDir,
                         const IVolumeMaterial* material) {
  constexpr bool isDense = Mode == SympyStepMode::Dense;

  double h = state.stepSize.value() * propDir;

  const double initialH = h;

  const Vector3 pos = stepper.position(state);
  const Vector3 dir = stepper.direction(state);
  const double t = stepper.time(state);
  const double qop = stepper.qOverP(state);
  const double pabs = stepper.absoluteMomentum(state);
  const double m = stepper.particleHypothesis(state).mass();

  if constexpr (isDense) {
    if (material != nullptr && pabs < state.options.dense.momentumCutOff) {
      return EigenStepperError::StepInvalid;
    }
  }

  const auto getB = [&](std::span<const double, 3> p) -> Result<Vector3> {
    return stepper.getField(state, {p[0], p[1], p[2]});
  };

  if (!state.field.has_value()) {
    auto fieldRes = stepper.getField(state, pos);
    if (!fieldRes.ok()) {
      return fieldRes.error();
    }
    state.field = *fieldRes;
  }
  // The kernel reads the start field before writing the last sample, so the
  // two can share storage; a copy would sit on the loop-carried chain.
  Vector3& lastField = *state.field;

  // Read once: the kernel writes through spans that point into `state`, so a
  // later read would be ordered behind all of its stores.
  const double stepTolerance = state.options.stepTolerance;

  const auto calcStepSizeScaling = [&](const double errorEstimate_) -> double {
    // For details about these values see ATL-SOFT-PUB-2009-001
    constexpr double lower = 0.25;
    constexpr double upper = 4.0;
    // This is given by the order of the Runge-Kutta method
    constexpr double exponent = 0.25;

    double x = stepTolerance / errorEstimate_;

    if constexpr (exponent == 0.25) {
      // This is 3x faster than std::pow
      x = std::sqrt(std::sqrt(x));
    } else {
      x = std::pow(x, exponent);
    }

    return std::clamp(x, lower, upper);
  };

  std::size_t nStepTrials = 0;
  double errorEstimate = 0.;

  while (true) {
    ++nStepTrials;
    ++state.statistics.nAttemptedSteps;

    // For details about the factor 4 see ATL-SOFT-PUB-2009-001
    std::error_code fieldError;
    const std::span<const double, 3> startPos(pos.data(), 3);
    const std::span<const double, 3> startDir(dir.data(), 3);
    const std::span<double, 3> endPos(
        state.pars.template segment<3>(eFreePos0).data(), 3);
    const std::span<double, 3> endDir(
        state.pars.template segment<3>(eFreeDir0).data(), 3);
    Rk4Status status{};
    if constexpr (isDense) {
      const std::span<double, 8> derivative(state.derivative.data(), 8);
      const std::span<double> jac =
          state.covTransport ? std::span<double>(state.jacToGlobal.data(),
                                                 state.jacToGlobal.size())
                             : std::span<double>();
      if (material == nullptr) {
        status = rk4_vacuum(startPos, startDir, t, h, qop, m, pabs,
                            std::span<const double, 3>(state.field->data(), 3),
                            getB, errorEstimate, 4 * stepTolerance, fieldError,
                            endPos, state.pars[eFreeTime], endDir,
                            std::span<double, 3>(lastField.data(), 3),
                            derivative, jac);
      } else {
        status = sympyDenseStep(stepper, state, *material, h, 4 * stepTolerance,
                                errorEstimate, lastField, fieldError, jac);
      }
    } else if constexpr (Mode == SympyStepMode::VacuumJac) {
      status =
          rk4_vacuum_jac(startPos, startDir, t, h, qop, m, pabs,
                         std::span<const double, 3>(state.field->data(), 3),
                         getB, errorEstimate, 4 * stepTolerance, fieldError,
                         endPos, state.pars[eFreeTime], endDir,
                         std::span<double, 3>(lastField.data(), 3),
                         std::span<double, 8>(state.derivative.data(), 8),
                         std::span<double>(state.jacToGlobal.data(),
                                           state.jacToGlobal.size()));
    } else {
      // No jacobian, so no path derivatives either.
      status =
          rk4_vacuum_nojac(startPos, startDir, t, h, qop, m, pabs,
                           std::span<const double, 3>(state.field->data(), 3),
                           getB, errorEstimate, 4 * stepTolerance, fieldError,
                           endPos, state.pars[eFreeTime], endDir,
                           std::span<double, 3>(lastField.data(), 3));
    }
    if (status == Rk4Status::FieldError) {
      return fieldError;
    }
    // Protect against division by zero
    errorEstimate = std::max(1e-20, errorEstimate);

    if (status == Rk4Status::Accepted) {
      break;
    }

    ++state.statistics.nRejectedSteps;

    const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
    h *= stepSizeScaling;

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(h) < std::abs(state.options.stepSizeCutOff)) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > state.options.maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
  }

  state.pathAccumulated += h;
  ++state.nSteps;
  state.nStepTrials += nStepTrials;

  ++state.statistics.nSuccessfulSteps;
  if (propDir != Direction::fromScalarZeroAsPositive(initialH)) {
    ++state.statistics.nReverseSteps;
  }
  state.statistics.pathLength += h;
  state.statistics.absolutePathLength += std::abs(h);

  const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
  const double nextAccuracy = std::abs(h * stepSizeScaling);
  const double previousAccuracy = std::abs(state.stepSize.accuracy());
  const double initialStepLength = std::abs(initialH);
  if (nextAccuracy < initialStepLength || nextAccuracy > previousAccuracy) {
    state.stepSize.setAccuracy(nextAccuracy);
  }

  if constexpr (isDense) {
    if (material != nullptr || !state.materialEffectsAccumulator.isVacuum()) {
      if (state.materialEffectsAccumulator.isVacuum()) {
        state.materialEffectsAccumulator.initialize(
            state.options.maxXOverX0Step, stepper.particleHypothesis(state),
            pabs);
      }

      Material mat =
          material != nullptr ? material->material(pos) : Material::Vacuum();

      state.materialEffectsAccumulator.accumulate(mat, propDir * h, qop,
                                                  stepper.qOverP(state));
    }
  }

  return h;
}

}  // namespace Acts::detail
