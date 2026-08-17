// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// The body of the vacuum step, shared by the two translation units that
// instantiate it -- one for covariance transport, one without.  Each sees
// exactly one instantiation and only the kernel that instantiation calls, so it
// gets its own stack frame and register allocation and clang inlines its
// kernel, which it does not do when both flavours sit in one translation unit.

#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>
#include <cmath>
#include <span>

// The kernel has to end up inlined into the step: the two want to be one basic
// block so that the scheduler can hide the jacobian transport in the value
// path's slack, and that is worth 8% with covariance transport.  Reached
// through a template instantiation clang stops doing it on its own, so pin the
// decision on.  GNU C++ allows an attribute between the return type and the
// declarator, which keeps this out of the generator.
#define ACTS_SYMPY_ALWAYS_INLINE(name) __attribute__((always_inline)) name
#define rk4_vacuum_jac ACTS_SYMPY_ALWAYS_INLINE(rk4_vacuum_jac)
#define rk4_vacuum_nojac ACTS_SYMPY_ALWAYS_INLINE(rk4_vacuum_nojac)
#include "codegen/sympy_stepper_math.hpp"
#undef rk4_vacuum_jac
#undef rk4_vacuum_nojac

namespace Acts::detail {

template <bool WithJac>
Result<double> sympyVacuumStep(const SympyStepper& stepper,
                               SympyStepper::State& state, Direction propDir) {
  double h = state.stepSize.value() * propDir;

  const double initialH = h;

  const Vector3 pos = stepper.position(state);
  const Vector3 dir = stepper.direction(state);
  const double t = stepper.time(state);
  const double qop = stepper.qOverP(state);
  const double pabs = stepper.absoluteMomentum(state);
  const double m = stepper.particleHypothesis(state).mass();

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
  // The kernel writes nothing until it has accepted a trial, so its last field
  // sample can go straight into the state instead of into a local that is
  // copied over afterwards.  That copy sat on the loop-carried chain: the field
  // at the end of one step is the first thing the next one multiplies.
  Vector3& lastField = *state.field;

  // `state.options.stepTolerance` is read once per trial to form the error
  // bound, and again afterwards to scale the step size for the next one. The
  // second read has to be ordered behind every store the kernel makes through
  // its output spans, because those point into `state` and the compiler cannot
  // prove they miss `state.options`. That puts the divide and the two square
  // roots of the feedback -- the longest dependency chain in the step, and one
  // that closes into the next step -- after the jacobian transport instead of
  // alongside it. Reading it once up front removes that ordering edge.
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

  // If a future path changes q/p without refreshing dt/ds, the times come out
  // silently wrong.  Catch it where it would happen rather than downstream.
  // Guarded on NDEBUG rather than left to `assert`: ACTS' assert_include shim
  // keeps asserts live in release builds, and this one costs the very sqrt and
  // divide the cached value exists to avoid.
#ifndef NDEBUG
  assert(std::abs(state.dtds - std::sqrt(1 + m * m / (pabs * pabs))) <
             1e-12 * state.dtds &&
         "cached dt/ds is stale: q/p changed without refreshing it");
#endif

  std::size_t nStepTrials = 0;
  double errorEstimate = 0.;

  while (true) {
    ++nStepTrials;
    ++state.statistics.nAttemptedSteps;

    // For details about the factor 4 see ATL-SOFT-PUB-2009-001
    bool accepted = false;
    std::error_code fieldError;
    const std::span<const double, 3> startPos(pos.data(), 3);
    const std::span<const double, 3> startDir(dir.data(), 3);
    const std::span<double, 3> endPos(
        state.pars.template segment<3>(eFreePos0).data(), 3);
    const std::span<double, 3> endDir(
        state.pars.template segment<3>(eFreeDir0).data(), 3);
    // A status code rather than a `Result<bool>`: a variant returned across
    // the kernel boundary is read back through the stack on the accepted
    // path, which is every path that matters.
    int status = 0;
    if constexpr (WithJac) {
      status = rk4_vacuum_jac(
          startPos, startDir, t, h, qop, m, pabs,
          std::span<const double, 3>(state.field->data(), 3), getB,
          errorEstimate, 4 * stepTolerance, fieldError, endPos,
          state.pars[eFreeTime], endDir,
          std::span<double, 3>(lastField.data(), 3), state.dtds,
          std::span<double, 8>(state.derivative.data(), 8),
          std::span<double>(state.jacToGlobal.data(), state.jacToGlobal.size()));
    } else {
      // The caller has already decided there is no jacobian, so neither the
      // jacobian span nor the path derivatives are formed at all.
      status = rk4_vacuum_nojac(
          startPos, startDir, t, h, qop, m, pabs,
          std::span<const double, 3>(state.field->data(), 3), getB,
          errorEstimate, 4 * stepTolerance, fieldError, endPos,
          state.pars[eFreeTime], endDir,
          std::span<double, 3>(lastField.data(), 3), state.dtds);
    }
    if (status == 2) {
      return fieldError;
    }
    accepted = status != 0;
    // Protect against division by zero
    errorEstimate = std::max(1e-20, errorEstimate);

    if (accepted) {
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

  return h;
}

}  // namespace Acts::detail
