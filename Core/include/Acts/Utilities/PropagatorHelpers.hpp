// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>

namespace Acts {
namespace PropagatorHelpers {

template <typename propagator_state_t, typename stepper_t>
inline auto charge(const propagator_state_t &state, const stepper_t &stepper) {
  assert(state.options.absCharge >= 0 && stepper.qop(state.stepping) != 0);
  return (stepper.qop(state.stepping) >= 0 ? +1 : -1) * state.options.absCharge;
}

template <typename propagator_state_t, typename stepper_t>
inline auto absoluteMomentum(const propagator_state_t &state,
                             const stepper_t &stepper) {
  auto q = charge(state, stepper);
  return (q == 0 ? 1 : q) / stepper.qop(state.stepping);
}

template <typename propagator_state_t, typename stepper_t>
inline Vector3 momentum(const propagator_state_t &state,
                        const stepper_t &stepper) {
  return absoluteMomentum(state, stepper) * stepper.direction(state.stepping);
}

}  // namespace PropagatorHelpers
}  // namespace Acts
