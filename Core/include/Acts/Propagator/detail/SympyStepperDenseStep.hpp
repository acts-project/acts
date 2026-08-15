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
#include "Acts/Propagator/detail/SympyStepperStatus.hpp"
#include "Acts/Utilities/Result.hpp"

#include <span>
#include <system_error>

namespace Acts {

class IVolumeMaterial;

namespace detail {

/// @brief A whole Runge-Kutta step of a track that is inside material
///
/// The whole dense path of `SympyStepper::step` -- momentum cut-off, kernel
/// branch and material accumulation -- so that none of it shares a stack frame
/// or a register allocation with the vacuum path.  Only called when
/// `state.options.doDense` and there is material in play.
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
/// @param [out] fieldErr the error of a failed field lookup
/// @param [in,out] jac the bound-to-free jacobian, empty to skip transport
///
/// @return whether the step was accepted, rejected or hit a field error
Rk4Status sympyDenseStep(const SympyStepper& stepper,
                         SympyStepper::State& state,
                         const IVolumeMaterial& material, double h,
                         double errTol, double& errorEstimate,
                         Vector3& lastField, std::error_code& fieldErr,
                         std::span<double> jac);

}  // namespace detail
}  // namespace Acts
