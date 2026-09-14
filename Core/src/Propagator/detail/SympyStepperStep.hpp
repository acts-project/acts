// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

class IVolumeMaterial;

namespace detail {

/// What a step has to carry: the dense path, or the vacuum one specialised on
/// covariance transport so the kernel does not test an empty jacobian span on
/// every trial.
enum class SympyStepMode {
  VacuumNoJac,
  VacuumJac,
  Dense,
};

/// @brief A whole Runge-Kutta step
///
/// One body for all three modes, instantiated once per mode in a translation
/// unit of its own; sharing one costs more than the branches it saves.
///
/// @param [in] stepper the stepper, for field access and accessors
/// @param [in,out] state the stepper state
/// @param [in] propDir the propagation direction
/// @param [in] material the volume material, `Dense` only, and null there when
///        just the accumulator still has material to flush
///
/// @return the step length taken, or an error
template <SympyStepMode Mode>
Result<double> sympyStep(const SympyStepper& stepper,
                         SympyStepper::State& state, Direction propDir,
                         const IVolumeMaterial* material);

extern template Result<double> sympyStep<SympyStepMode::VacuumNoJac>(
    const SympyStepper&, SympyStepper::State&, Direction,
    const IVolumeMaterial*);
extern template Result<double> sympyStep<SympyStepMode::VacuumJac>(
    const SympyStepper&, SympyStepper::State&, Direction,
    const IVolumeMaterial*);
extern template Result<double> sympyStep<SympyStepMode::Dense>(
    const SympyStepper&, SympyStepper::State&, Direction,
    const IVolumeMaterial*);

}  // namespace detail
}  // namespace Acts
