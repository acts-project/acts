// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/stepping_config.hpp"
#include "detray/tracks/ray.hpp"

// System include(s)
#include <iomanip>
#include <sstream>
#include <string>

namespace detray::stepping {

/// Print basic information about the state of a stepping stream @param state
template <typename state_type, concepts::scalar scalar_t>
DETRAY_HOST inline std::string print_state(const state_type &state,
                                           const scalar_t dist) {
  // Gathers stepping information across stepper update calls
  std::stringstream debug_stream{};
  // Column width in output
  constexpr int cw{20};

  debug_stream << std::left << std::setw(cw) << "Step size" << state.step_size()
               << std::endl;
  debug_stream << std::setw(cw) << "Dist to sf." << dist << std::endl;
  debug_stream << std::setw(cw) << "Path length" << state.path_length()
               << std::endl;

  switch (state.direction()) {
    using enum step::direction;
    case e_forward:
      debug_stream << std::setw(cw) << "direction" << "forward" << std::endl;
      break;
    case e_unknown:
      debug_stream << std::setw(cw) << "direction" << "unknown" << std::endl;
      break;
    case e_backward:
      debug_stream << std::setw(cw) << "direction" << "backward" << std::endl;
      break;
    default:
      break;
  }

  auto pos = state().pos();

  debug_stream << std::setw(cw)
               << "Pos: " << "[r = " << math::hypot(pos[0], pos[1])
               << ", z = " << pos[2] << "]" << std::endl;
  debug_stream << std::setw(cw) << "Tangent:" << detail::ray(state());

  return debug_stream.str();
}

/// Print basic information about the state of a stepping stream @param state
template <typename state_type, concepts::scalar scalar_t>
DETRAY_HOST inline std::string print_state(const state_type &state,
                                           const std::size_t n_trials,
                                           const scalar_t step_scalor) {
  // Gathers stepping information across stepper update calls
  std::stringstream debug_stream{};
  // Column width in output
  constexpr int cw{20};

  // Remove trailing newlines
  debug_stream << std::left << std::setw(cw)
               << "Step size: " << state.step_size() << std::endl;
  debug_stream << std::setw(cw) << "no. RK adjustments: " << n_trials
               << std::endl;
  debug_stream << std::setw(cw) << "Step size scale factor: " << step_scalor;

  return debug_stream.str();
}

}  // namespace detray::stepping
