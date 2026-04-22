// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/navigation_config.hpp"

// System include(s)
#include <iomanip>
#include <sstream>
#include <string>

namespace detray::navigation {

/// Print basic information about the state of a navigation stream @param state
template <typename state_type>
DETRAY_HOST inline std::string print_state(const state_type &state) {
  using detector_t = typename state_type::detector_type;
  using scalar_t = typename detector_t::scalar_type;

  // Gathers navigation information across navigator update calls
  std::stringstream debug_stream{};
  // Column width in output
  constexpr int cw{20};

  debug_stream << std::left << std::setw(cw) << "Volume:" << state.volume()
               << std::endl;

  debug_stream << std::setw(cw) << std::boolalpha
               << "is alive:" << state.is_alive() << std::endl;
  debug_stream << std::noboolalpha;

  // Navigation direction
  debug_stream << std::setw(cw) << "direction:";
  debug_stream << state.direction() << std::endl;

  // Navigation status
  debug_stream << std::setw(cw) << "status:";
  debug_stream << state.status() << std::endl;

  // Navigation trust level
  debug_stream << std::setw(cw) << "trust:";
  debug_stream << state.trust_level() << std::endl;

  // Number of reachable candidates
  debug_stream << std::setw(cw) << "No. reachable:" << state.n_candidates()
               << std::endl;

  // Current surface
  debug_stream << std::setw(cw) << "current object:";
  if (state.is_on_surface()) {
    // If "exit" is called twice, the state has been cleared
    debug_stream << state.geometry_identifier() << std::endl;
  } else if (state.status() == status::e_exit) {
    debug_stream << "exited" << std::endl;
  } else {
    debug_stream << "undefined" << std::endl;
  }

  // Next surface
  if (!state.candidates().empty()) {
    debug_stream << std::setw(cw) << "next object:";
    if (state.n_candidates() == 0u) {
      debug_stream << "exhausted" << std::endl;
    } else {
      debug_stream << state.next_surface().identifier() << std::endl;
    }
  }

  // Distance to next
  debug_stream << std::setw(cw) << "distance to next:";
  if (!state.cache_exhausted() && state.is_on_surface()) {
    debug_stream << "on obj (within tol)" << std::endl;
  } else if (state.cache_exhausted()) {
    debug_stream << "no target" << std::endl;
  } else {
    debug_stream << state() / detray::unit<scalar_t>::mm << " mm" << std::endl;
  }

  // Current external mask tolerance
  debug_stream << std::setw(cw) << "ext. mask tol.:"
               << state.external_tol() / detray::unit<scalar_t>::mm << " mm"
               << std::endl;

  return debug_stream.str();
}

/// Print candidate and configuration information of a navigation state
///
/// @param state the state object of the navigation stream
/// @param cfg the navigation configuration object
/// @param track_pos the current track position
/// @param track_dir the current track direction
template <typename state_type, concepts::point3D point3_t,
          concepts::vector3D vector3_t>
DETRAY_HOST inline std::string print_candidates(const state_type &state,
                                                const navigation::config &cfg,
                                                const point3_t &track_pos,
                                                const vector3_t &track_dir) {
  using detector_t = typename state_type::detector_type;
  using geo_ctx_t = typename detector_t::geometry_context;
  using scalar_t = typename detector_t::scalar_type;

  // Gathers navigation information across navigator update calls
  std::stringstream debug_stream{};
  // Column width in output
  constexpr int cw{20};

  debug_stream << std::left << std::setw(cw) << "Overstep tol.:"
               << cfg.intersection.overstep_tolerance /
                      detray::unit<scalar_t>::um
               << " um" << std::endl;

  debug_stream << std::setw(cw) << "Track:"
               << "pos: [r = " << vector::perp(track_pos)
               << ", z = " << track_pos[2] << "]," << std::endl;

  debug_stream << std::setw(cw) << " "
               << "dir: [" << track_dir[0] << ", " << track_dir[1] << ", "
               << track_dir[2] << "]" << std::endl;

  debug_stream << "Surface candidates: " << std::endl;

  for (const auto &sf_cand : state) {
    debug_stream << std::left << std::setw(6) << "-> " << sf_cand;

    assert(!sf_cand.surface().identifier().is_invalid());

    // Use additional debug information that was gathered on the cand.
    if constexpr (state_type::value_type::contains_pos()) {
      const auto &local = sf_cand.local();
      if (!sf_cand.surface().identifier().is_invalid()) {
        point3_t pos = geometry::surface{state.detector(), sf_cand.surface()}
                           .local_to_global(geo_ctx_t{}, local, track_dir);
        debug_stream << " glob: [r = " << vector::perp(pos)
                     << ", z = " << pos[2] << "]" << std::endl;
      } else {
        debug_stream << "Invalid identifier" << std::endl;
      }
    } else {
      debug_stream << std::endl;
    }
  }

  return debug_stream.str();
}

}  // namespace detray::navigation
