// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/inspectors.hpp"
#include "detray/test/utils/step_tracer.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Covfie include(s)
#include <covfie/core/field.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s)
#include <stdexcept>

namespace detray {

// These types are identical in host and device code for all bfield types
using test_algebra = test::algebra;
using scalar = dscalar<test_algebra>;
using vector3 = dvector3D<test_algebra>;
using point3 = dpoint3D<test_algebra>;
using test_track = free_track_parameters<test_algebra>;
using free_matrix_t = free_matrix<test_algebra>;

constexpr std::size_t cache_size{navigation::default_cache_size};

// Navigator
template <typename detector_t, typename inspector_t>
using navigator_w_insp_t =
    caching_navigator<detector_t, cache_size, inspector_t>;
template <typename detector_t>
using navigator_t = navigator_w_insp_t<detector_t, navigation::void_inspector>;
template <typename detector_t>
using intersection_t = typename navigator_t<detector_t>::intersection_type;

// Stepper
using constraints_t = constrained_step<scalar>;
template <typename bfield_view_t>
using rk_stepper_t = rk_stepper<bfield_view_t, test_algebra, constraints_t>;

// Track generator
using generator_t = uniform_track_generator<test_track>;

/// Test tolerance
constexpr scalar is_close{1e-4f};

/// Test configuration
struct propagator_test_config {
  generator_t::configuration track_generator;
  propagation::config propagation;
};

// Assemble actor chain type
using step_tracer_host_t = actor::step_tracer<test_algebra, vecmem::vector>;
using step_tracer_device_t =
    actor::step_tracer<test_algebra, vecmem::device_vector>;
using pathlimit_aborter_t = actor::pathlimit_aborter<scalar>;
using parameter_setter_t = actor::parameter_setter<test_algebra>;
using actor_chain_host_t = actor_chain<
    pathlimit_aborter_t,
    actor::parameter_updater<test_algebra,
                             actor::pointwise_material_interactor<test_algebra>,
                             step_tracer_host_t>>;
using actor_chain_device_t = actor_chain<
    pathlimit_aborter_t,
    actor::parameter_updater<test_algebra,
                             actor::pointwise_material_interactor<test_algebra>,
                             step_tracer_device_t>>;

/// Precompute the tracks
template <typename track_generator_t = uniform_track_generator<test_track>>
inline auto generate_tracks(
    vecmem::memory_resource *mr,
    const typename track_generator_t::configuration &cfg = {}) {
  // Track collection
  vecmem::vector<typename track_generator_t::track_type> tracks(mr);

  // Iterate through uniformly distributed momentum directions
  for (auto track : track_generator_t{cfg}) {
    // Put it into vector of trajectories
    tracks.push_back(track);
  }

  return tracks;
}

/// Test function for propagator on the host
template <typename bfield_bknd_t, typename host_detector_t>
inline auto run_propagation_host(vecmem::memory_resource *mr,
                                 const host_detector_t &det,
                                 const propagation::config &cfg,
                                 const covfie::field<bfield_bknd_t> &field,
                                 const vecmem::vector<test_track> &tracks)
    -> vecmem::jagged_vector<step_record<test_algebra>> {
  // Construct propagator from stepper and navigator
  using host_stepper_t =
      rk_stepper_t<typename covfie::field<bfield_bknd_t>::view_t>;
  using host_navigator_t =
      navigator_w_insp_t<host_detector_t, navigation::print_inspector>;
  using propagator_host_t =
      propagator<host_stepper_t, host_navigator_t, actor_chain_host_t>;

  propagator_host_t p{cfg};

  // Create vector for track recording
  vecmem::jagged_vector<step_record<test_algebra>> host_steps(mr);

  for (const auto &trk : tracks) {
    // Create the propagator state
    step_tracer_host_t::state tracer_state{*mr};
    typename pathlimit_aborter_t::state pathlimit_state{
        cfg.stepping.path_limit};
    actor::parameter_updater_state<test_algebra> updater_state{cfg};
    actor::pointwise_material_interactor<test_algebra>::state
        interactor_state{};

    auto actor_states = detray::tie(tracer_state, pathlimit_state,
                                    updater_state, interactor_state);

    typename covfie::field<bfield_bknd_t>::view_t bfield_view(field);
    typename propagator_host_t::state state(trk, bfield_view, det);

    state.stepping().template set_constraint<step::constraint::e_accuracy>(
        cfg.stepping.step_constraint);

    // Run propagation
    if (!p.propagate(state, actor_states)) {
      DETRAY_FATAL_HOST(state.navigation().inspector().to_string());
      throw std::runtime_error("Host propagation failed");
    } else if (tracer_state.get_step_data().empty()) {
      DETRAY_FATAL_HOST(state.navigation().inspector().to_string());
      throw std::runtime_error(
          "Host propagation did not record reference data correctly");
    }

    host_steps.push_back(std::move(tracer_state.get_step_data()));
  }

  return host_steps;
}

/// Compare the results between host and device propagation
inline void compare_propagation_results(
    const vecmem::jagged_vector<step_record<test_algebra>> &host_steps,
    const vecmem::jagged_vector<step_record<test_algebra>> &device_steps) {
  // Make sure the same number of tracks were tested on both backends
  ASSERT_EQ(host_steps.size(), device_steps.size());

  // Compare the positions and pathlengths
  for (unsigned int i = 0u; i < host_steps.size(); i++) {
    EXPECT_TRUE(host_steps[i].size() > 0u)
        << "No surfaces forund on track " << i;

    // Recorded as many positions as path lengths
    EXPECT_EQ(host_steps[i].size(), device_steps[i].size());

    for (unsigned int j = 0u;
         j < math::min(host_steps[i].size(), device_steps[i].size()); j++) {
      // Compare recorded path lengths along track
      const auto &host_step = host_steps[i][j];
      const auto &device_step = device_steps[i][j];

      EXPECT_NEAR(host_step.path_length, device_step.path_length,
                  host_step.path_length * is_close)
          << "ERROR: Path length at track " << i << " step " << j << std::endl;

      // Compare recorded positions along track
      const point3 &host_pos = host_step.free_params.pos();
      const point3 &device_pos = device_step.free_params.pos();

      auto relative_error = static_cast<point3>(1.f / host_step.path_length *
                                                (host_pos - device_pos));

      EXPECT_NEAR(vector::norm(relative_error), 0.f, is_close)
          << "ERROR: Position at track " << i << " step " << j << ": ["
          << host_pos[0] << ", " << host_pos[1] << ", " << host_pos[2]
          << "] (host), [" << device_pos[0] << ", " << device_pos[1] << ", "
          << device_pos[2] << "] (device)" << std::endl;

      // Compare the Jacobians
      const free_matrix_t &host_J = host_step.jacobian;
      const free_matrix_t &device_J = device_step.jacobian;

      for (std::size_t row = 0u; row < e_free_size; row++) {
        for (std::size_t col = 0u; col < e_free_size; col++) {
          scalar host_val = getter::element(host_J, row, col);

          scalar device_val = getter::element(device_J, row, col);

          ASSERT_NEAR((host_val - device_val) / host_step.path_length, 0.f,
                      is_close)
              << "ERROR: matrix element mismatch at row " << row << ", col "
              << col << std::endl;
        }
      }
    }
  }
}

}  // namespace detray
