// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Utilities/Result.hpp"

#include <span>

namespace Acts {

class IVolumeMaterial;

namespace detail {

/// @brief A whole Runge-Kutta step of a track that is inside material
///
/// Specialised on covariance transport, each specialisation in a translation
/// unit of its own, for the same reason as the vacuum step.
///
/// The dense path of `SympyStepper::step`, moved out of it wholesale.  Keeping
/// only the kernel in its own translation unit was not enough: the momentum
/// cut-off, the branch and the material accumulation that stayed behind still
/// shared a stack frame and a register allocation with the vacuum path, which
/// is the one that runs almost everywhere.
///
/// @param [in] stepper the stepper, for field access and accessors
/// @param [in,out] state the stepper state
/// @param [in] propDir the propagation direction
/// @param [in] material the volume material, may be null when only the
///        accumulator still has material to flush
///
/// @return the step length taken, or an error
template <bool WithJac>
Result<double> sympyDenseStepFull(const SympyStepper& stepper,
                                  SympyStepper::State& state, Direction propDir,
                                  const IVolumeMaterial* material);

}  // namespace detail
}  // namespace Acts
