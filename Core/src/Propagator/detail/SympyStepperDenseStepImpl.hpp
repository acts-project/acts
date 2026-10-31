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
#include "Acts/Material/detail/EnergyLossConstants.hpp"
#include "Acts/Propagator/detail/SympyStepperDenseStep.hpp"

#include <cmath>

#include "codegen/sympy_stepper_math.hpp"

namespace Acts {

template <bool WithJac>
detail::Rk4Status detail::sympyDenseStep(
    const SympyStepper& stepper, SympyStepper::State& state,
    const IVolumeMaterial& material, double h, double errTol,
    double& errorEstimate, Vector3& lastField, std::error_code& fieldErr,
    std::span<double> jac) {
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

  if constexpr (WithJac) {
    const auto getGGradient = [&](std::span<const double, 3> p,
                                  double l) -> double {
      if (!state.options.dense.includeGradient ||
          stepper.particleHypothesis(state).extractMomentum(l) <
              state.options.dense.momentumCutOff) {
        return 0.;
      }
      const MaterialSlab slab(material.material({p[0], p[1], p[2]}),
                              1.0f * UnitConstants::mm);
      if (state.options.dense.meanEnergyLoss) {
        return deriveEnergyLossMeanQOverP(slab, absPdg, static_cast<float>(m),
                                          static_cast<float>(l), absQ);
      }
      return deriveEnergyLossModeQOverP(slab, absPdg, static_cast<float>(m),
                                        static_cast<float>(l), absQ);
    };
    const auto status = rk4_dense_jac(
        std::span<const double, 3>(pos.data(), 3),
        std::span<const double, 3>(dir.data(), 3), t, h, qop, m, q, pabs,
        std::span<const double, 3>(state.field->data(), 3), getB, getG,
        getGGradient, errorEstimate, errTol, fieldErr,
        std::span<double, 3>(state.pars.segment<3>(eFreePos0).data(), 3),
        state.pars[eFreeTime],
        std::span<double, 3>(state.pars.segment<3>(eFreeDir0).data(), 3),
        state.pars[eFreeQOverP], std::span<double, 3>(lastField.data(), 3),
        std::span<double, 8>(state.derivative.data(), 8), jac);

    if (status == Rk4Status::Accepted && state.options.dense.includeGradient) {
      // A smooth variational equation misses the change in crossing time at
      // a discontinuity of the stopping power. For a boundary q/p = constant,
      // the saltation matrix is the identity except for its q/p diagonal:
      // f_after / f_before = (dE/ds)_after / (dE/ds)_before.
      // Apply it at the end of the accepted crossing step. The RK error
      // controller resolves the jump; the placement error vanishes with that
      // step's length. No correction is needed on rejected trials.
      const double newPabs = stepper.absoluteMomentum(state);
      const auto correctAt = [&](double threshold) {
        if ((pabs < threshold) == (newPabs < threshold)) {
          return;
        }
        const double boundaryQop = q / threshold;
        // Loss functions use float internally. Sample far enough from the
        // boundary to select each branch reliably, then extrapolate its
        // smooth derivative back to the boundary.
        constexpr double relativeOffset = 1e-5;
        const auto boundaryLoss = [&](double offset) {
          const double sample = boundaryQop * (1 + offset);
          const std::span<const double, 3> endPos(
              state.pars.template segment<3>(eFreePos0).data(), 3);
          return getG(endPos, sample) +
                 getGGradient(endPos, sample) * (boundaryQop - sample);
        };
        const double lowMomentumLoss = boundaryLoss(relativeOffset);
        const double highMomentumLoss = boundaryLoss(-relativeOffset);
        // Vacuum, or a user momentum cutoff, has no regular crossing flow.
        if (lowMomentumLoss <= 0 || highMomentumLoss <= 0) {
          return;
        }
        const double factor = newPabs < pabs
                                  ? lowMomentumLoss / highMomentumLoss
                                  : highMomentumLoss / lowMomentumLoss;
        Eigen::Map<BoundToFreeMatrix> matrix(jac.data());
        matrix.row(eFreeQOverP) *= factor;
        // Preserve the other rows in the plain Jacobian while restoring the
        // kernel's scaling by the current q/p and its transported derivative.
        matrix.template block<eFreeQOverP, 1>(0, eBoundQOverP) /= factor;
      };
      correctAt(BetheDensityEffectBetaGammaThreshold * m);
      if (absPdg == PdgParticle::eMuon) {
        for (double energy : {MuonRadiativeLowEnergyThreshold,
                              MuonRadiativeHighEnergyThreshold}) {
          if (energy > m) {
            correctAt(std::sqrt(energy * energy - m * m));
          }
        }
      }
    }
    return status;
  } else {
    return rk4_dense_nojac(
        std::span<const double, 3>(pos.data(), 3),
        std::span<const double, 3>(dir.data(), 3), t, h, qop, m, q, pabs,
        std::span<const double, 3>(state.field->data(), 3), getB, getG,
        errorEstimate, errTol, fieldErr,
        std::span<double, 3>(state.pars.segment<3>(eFreePos0).data(), 3),
        state.pars[eFreeTime],
        std::span<double, 3>(state.pars.segment<3>(eFreeDir0).data(), 3),
        state.pars[eFreeQOverP], std::span<double, 3>(lastField.data(), 3));
  }
}

}  // namespace Acts
