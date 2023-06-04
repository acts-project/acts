// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

namespace Acts {
namespace PropagatorHelpers {

template <typename propagator_state_t, typename stepper_t>
inline auto charge(const propagator_state_t &state, const stepper_t &stepper) {
  return std::copysign(stepper.qop(state.stepping), state.options.absCharge);
}

template <typename propagator_state_t, typename stepper_t>
inline auto absoluteMomentum(const propagator_state_t &state,
                             const stepper_t &stepper) {
  return state.options.absCharge / std::abs(stepper.qop(state.stepping));
}

template <typename propagator_state_t, typename stepper_t>
inline auto momentum(const propagator_state_t &state,
                     const stepper_t &stepper) {
  return absoluteMomentum(state, stepper) * stepper.direction(state);
}

}  // namespace PropagatorHelpers
}  // namespace Acts
