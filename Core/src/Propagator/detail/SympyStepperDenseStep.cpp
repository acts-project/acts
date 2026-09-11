// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/SympyStepperDenseStep.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include <cmath>

#include "codegen/sympy_stepper_math.hpp"

namespace Acts::detail {

Rk4Status sympyDenseStep(const SympyStepper& stepper,
                         SympyStepper::State& state,
                         const IVolumeMaterial& material, double h,
                         double errTol, double& errorEstimate,
                         Vector3& lastField, std::error_code& fieldErr,
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

  return rk4_dense(
      std::span<const double, 3>(pos.data(), 3),
      std::span<const double, 3>(dir.data(), 3), t, h, qop, m, q, pabs,
      std::span<const double, 3>(state.field->data(), 3), getB, getG,
      errorEstimate, errTol, fieldErr,
      std::span<double, 3>(state.pars.segment<3>(eFreePos0).data(), 3),
      state.pars[eFreeTime],
      std::span<double, 3>(state.pars.segment<3>(eFreeDir0).data(), 3),
      state.pars[eFreeQOverP], std::span<double, 3>(lastField.data(), 3),
      std::span<double, 8>(state.derivative.data(), 8), jac);
}

}  // namespace Acts::detail
