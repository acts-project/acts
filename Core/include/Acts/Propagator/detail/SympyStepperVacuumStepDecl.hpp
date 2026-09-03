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

/// @brief A whole Runge-Kutta step of a track in vacuum
///
/// Specialised on covariance transport: the kernel is generated twice, with and
/// without the jacobian transport, and `SympyStepper::step` picks by
/// `state.covTransport`.  One kernel testing an empty jacobian span instead
/// splits itself into two basic blocks, and everything the jacobian half needs
/// is live across that edge -- 52 spill/reloads and 144 bytes of stack frame
/// that buy nothing, since the floating-point work is identical either way.
///
/// Each specialisation is instantiated in a translation unit of its own so that
/// it gets its own stack frame and register allocation and inlines its own
/// kernel; both in one unit costs more than the branch it saves.
///
/// @param [in] stepper the stepper, for field access and accessors
/// @param [in,out] state the stepper state
/// @param [in] propDir the propagation direction
/// @return the step length taken, or an error
template <bool WithJac>
Result<double> sympyVacuumStep(const SympyStepper& stepper,
                               SympyStepper::State& state, Direction propDir);

}  // namespace detail
}  // namespace Acts
