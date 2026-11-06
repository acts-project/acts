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

    // The stopping power jumps at the thresholds of the loss model. Crossing
    // one scales the q/p row of the jacobian by the ratio of the slopes after
    // and before it (saltation). This applies the factor at the step end, so
    // the error in the other rows is first order in the step length.
    if (status == Rk4Status::Accepted && state.options.dense.includeGradient) {
      const double newPabs = stepper.absoluteMomentum(state);
      const std::span<const double, 3> newPos(
          state.pars.segment<3>(eFreePos0).data(), 3);
      const auto saltate = [&](double threshold) {
        if ((pabs < threshold) == (newPabs < threshold)) {
          return;
        }
        // The loss is computed in float. Sample each side away from the
        // threshold and extrapolate linearly back to it.
        const double thresholdQop = q / threshold;
        const auto lossAt = [&](double offset) {
          const double l = thresholdQop * (1 + offset);
          return getG(newPos, l) + getGGradient(newPos, l) * (thresholdQop - l);
        };
        // A positive offset is on the low momentum side.
        const double offset = pabs < threshold ? 1e-5 : -1e-5;
        const double before = lossAt(offset);
        const double after = lossAt(-offset);
        // Vacuum or momentum cutoff
        if (before <= 0 || after <= 0) {
          return;
        }
        const double factor = after / before;
        Eigen::Map<BoundToFreeMatrix> J(jac.data());
        J.row(eFreeQOverP) *= factor;
        // Rows above q/p store J(i, q/p) * q/p / J(q/p, q/p), see
        // `toScaledBoundToFree`. Keep J(i, q/p) unchanged.
        J.block<eFreeQOverP, 1>(0, eBoundQOverP) /= factor;
      };
      saltate(BetheDensityEffectBetaGammaThreshold * m);
      if (absPdg == PdgParticle::eMuon) {
        for (double energy : {MuonRadiativeLowEnergyThreshold,
                              MuonRadiativeHighEnergyThreshold}) {
          saltate(std::sqrt(energy * energy - m * m));
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
