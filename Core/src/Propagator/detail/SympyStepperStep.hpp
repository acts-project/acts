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

/// Which path a step takes
enum class SympyStepMode {
  Vacuum,
  Dense,
};

/// @brief A whole Runge-Kutta step
///
/// One body for every combination, instantiated once per combination in a
/// translation unit of its own; sharing one costs more than the branches it
/// saves.
///
/// @tparam Mode which path the step takes
/// @tparam WithJac whether the jacobian is transported; specialising on it
///         keeps the kernel from testing an empty jacobian span on every trial
///
/// @param [in] stepper the stepper, for field access and accessors
/// @param [in,out] state the stepper state
/// @param [in] propDir the propagation direction
/// @param [in] material the volume material, @c Dense only, and null there
///        when just the accumulator still has material to flush
///
/// @return the step length taken, or an error
template <SympyStepMode Mode, bool WithJac>
Result<double> sympyStep(const SympyStepper& stepper,
                         SympyStepper::State& state, Direction propDir,
                         const IVolumeMaterial* material);

extern template Result<double> sympyStep<SympyStepMode::Vacuum, false>(
    const SympyStepper&, SympyStepper::State&, Direction,
    const IVolumeMaterial*);
extern template Result<double> sympyStep<SympyStepMode::Vacuum, true>(
    const SympyStepper&, SympyStepper::State&, Direction,
    const IVolumeMaterial*);
extern template Result<double> sympyStep<SympyStepMode::Dense, false>(
    const SympyStepper&, SympyStepper::State&, Direction,
    const IVolumeMaterial*);
extern template Result<double> sympyStep<SympyStepMode::Dense, true>(
    const SympyStepper&, SympyStepper::State&, Direction,
    const IVolumeMaterial*);

}  // namespace detail
}  // namespace Acts
