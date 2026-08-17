// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <span>

// One dense kernel per translation unit, inlined, for the same reason as the
// vacuum step: two of them in one function share a stack frame and a register
// allocation.
#define ACTS_SYMPY_ALWAYS_INLINE_D(name) __attribute__((always_inline)) name
#define rk4_dense ACTS_SYMPY_ALWAYS_INLINE_D(rk4_dense)
#define rk4_dense_jac ACTS_SYMPY_ALWAYS_INLINE_D(rk4_dense_jac)
#define rk4_dense_nojac ACTS_SYMPY_ALWAYS_INLINE_D(rk4_dense_nojac)
#include "codegen/sympy_stepper_math.hpp"
#undef rk4_dense
#undef rk4_dense_jac
#undef rk4_dense_nojac

namespace Acts::detail {
template <bool WithJac>
Result<bool> sympyDenseStep(const SympyStepper& stepper,
                            SympyStepper::State& state,
                            const IVolumeMaterial& material, double h,
                            double errTol, double& errorEstimate,
                            Vector3& lastField, std::span<double> jac) {
  const Vector3 pos = stepper.position(state);
  const Vector3 dir = stepper.direction(state);
  const double t = stepper.time(state);
  const double qop = stepper.qOverP(state);
  const double pabs = stepper.absoluteMomentum(state);
  const double m = stepper.particleHypothesis(state).mass();
  const double q = stepper.charge(state);
  const auto absQ = static_cast<float>(std::abs(q));
  const PdgParticle absPdg = stepper.particleHypothesis(state).absolutePdg();

  const auto getB = [&](std::span<const double, 3> p) {
    return stepper.getField(state, {p[0], p[1], p[2]});
  };

  const auto getG = [&](std::span<const double, 3> p, double l) -> double {
    if (const double newPabs =
            stepper.particleHypothesis(state).extractMomentum(l);
        newPabs < state.options.dense.momentumCutOff) {
      return 0.;
    }

    const MaterialSlab slab(material.material({p[0], p[1], p[2]}),
                            1.0f * UnitConstants::mm);
    // Unsigned: the step's own sign gives the energy back on a backward step,
    // matching `PointwiseMaterialInteraction`.
    if (state.options.dense.meanEnergyLoss) {
      return computeEnergyLossMean(slab, absPdg, static_cast<float>(m),
                                   static_cast<float>(l), absQ);
    }
    return computeEnergyLossMode(slab, absPdg, static_cast<float>(m),
                                 static_cast<float>(l), absQ);
  };

  // The caller already knows whether the jacobian is wanted, so the kernel
  // does not have to test an empty span on every trial.
  if constexpr (WithJac) {
    return rk4_dense_jac(
        std::span<const double, 3>(pos.data(), 3),
        std::span<const double, 3>(dir.data(), 3), t, h, qop, m, q, pabs,
        std::span<const double, 3>(state.field->data(), 3), getB, getG,
        errorEstimate, errTol,
        std::span<double, 3>(state.pars.segment<3>(eFreePos0).data(), 3),
        state.pars[eFreeTime],
        std::span<double, 3>(state.pars.segment<3>(eFreeDir0).data(), 3),
        state.pars[eFreeQOverP], std::span<double, 3>(lastField.data(), 3),
        std::span<double, 8>(state.derivative.data(), 8), jac);
  } else {
    return rk4_dense_nojac(
        std::span<const double, 3>(pos.data(), 3),
        std::span<const double, 3>(dir.data(), 3), t, h, qop, m, q, pabs,
        std::span<const double, 3>(state.field->data(), 3), getB, getG,
        errorEstimate, errTol,
        std::span<double, 3>(state.pars.segment<3>(eFreePos0).data(), 3),
        state.pars[eFreeTime],
        std::span<double, 3>(state.pars.segment<3>(eFreeDir0).data(), 3),
        state.pars[eFreeQOverP], std::span<double, 3>(lastField.data(), 3));
  }
}

template <bool WithJac>
Result<double> sympyDenseStepFull(const SympyStepper& stepper,
                                  SympyStepper::State& state, Direction propDir,
                                  const IVolumeMaterial* material) {
  double h = state.stepSize.value() * propDir;

  const double initialH = h;

  const Vector3 pos = stepper.position(state);
  const Vector3 dir = stepper.direction(state);
  const double t = stepper.time(state);
  const double qop = stepper.qOverP(state);
  const double pabs = stepper.absoluteMomentum(state);
  const double m = stepper.particleHypothesis(state).mass();

  if (state.options.doDense && material != nullptr &&
      pabs < state.options.dense.momentumCutOff) {
    return EigenStepperError::StepInvalid;
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
  // The kernels write nothing until they have accepted a trial, so the last
  // field sample can go straight into the state instead of into a local that
  // is copied over after the retry loop.  That copy sat on the loop-carried
  // chain: the field at the end of one step is the first thing the next one
  // multiplies.
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
    const std::span<double, 8> derivative(state.derivative.data(), 8);
    const std::span<double> jac =
        state.covTransport ? std::span<double>(state.jacToGlobal.data(),
                                               state.jacToGlobal.size())
                           : std::span<double>();
    if (!state.options.doDense || material == nullptr) {
      // A status code rather than a `Result<bool>`: a variant returned across
      // the kernel boundary is read back through the stack on the accepted
      // path, which is every path that matters.
      const int status = rk4_vacuum(
          startPos, startDir, t, h, qop, m, pabs,
          std::span<const double, 3>(state.field->data(), 3), getB,
          errorEstimate, 4 * stepTolerance, fieldError, endPos,
          state.pars[eFreeTime], endDir,
          std::span<double, 3>(lastField.data(), 3), derivative, jac);
      if (status == 2) {
        return fieldError;
      }
      accepted = status != 0;
    } else {
      Result<bool> res =
          sympyDenseStep<WithJac>(stepper, state, *material, h,
                                  4 * stepTolerance, errorEstimate, lastField,
                                  jac);
      if (!res.ok()) {
        return res.error();
      }
      accepted = *res;
    }
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

  // a dense step is the one place q/p moves, so the vacuum step's dt/ds has to
  // be brought back in line here
  {
    const double mDt = stepper.particleHypothesis(state).mass();
    const double pDt = stepper.absoluteMomentum(state);
    state.dtds = std::sqrt(1 + mDt * mDt / (pDt * pDt));
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

  if (state.options.doDense &&
      (material != nullptr || !state.materialEffectsAccumulator.isVacuum())) {
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

  return h;
}

}  // namespace Acts::detail
