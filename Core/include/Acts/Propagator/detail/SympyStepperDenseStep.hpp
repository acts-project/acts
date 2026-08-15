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
Result<double> sympyDenseStepFull(const SympyStepper& stepper,
                                  SympyStepper::State& state, Direction propDir,
                                  const IVolumeMaterial* material);

/// @brief One Runge-Kutta step through material
///
/// Kept in its own translation unit so that the dense kernel is not
/// instantiated next to SympyStepper::step, where it would share a stack frame
/// and a register allocation with the vacuum path.
///
/// @param [in] stepper the stepper, for field access
/// @param [in,out] state the stepper state, read for the start parameters and
///        written with the end parameters on success
/// @param [in] material the volume material
/// @param [in] h the step size to attempt
/// @param [in] errTol the tolerated error estimate
/// @param [out] errorEstimate the error estimate of the attempted step
/// @param [out] lastField the field at the last sampled point, to seed the
///        next step
/// @param [in,out] jac the bound-to-free jacobian, empty to skip transport
///
/// @return whether the step was accepted, or an error
Result<bool> sympyDenseStep(const SympyStepper& stepper,
                            SympyStepper::State& state,
                            const IVolumeMaterial& material, double h,
                            double errTol, double& errorEstimate,
                            Vector3& lastField, std::span<double> jac);

}  // namespace detail
}  // namespace Acts
