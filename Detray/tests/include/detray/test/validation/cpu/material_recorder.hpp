// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/navigation/detail/print_state.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/tracks/free_track_parameters.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/type_list.hpp"

// Detray IO include(s)
#include "detray/io/utils/create_path.hpp"
#include "detray/io/utils/file_handle.hpp"

// Test include(s)
#include "detray/test/validation/material_recording.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <array>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>

namespace detray::validation::host {

/// Run the propagation and record test data along the way
template <concepts::detector detector_t>
inline auto record_material(
    vecmem::memory_resource *host_mr, const detector_t &det,
    const propagation::config &cfg,
    const free_track_parameters<typename detector_t::algebra_type> &track) {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;

  using stepper_t = line_stepper<algebra_t>;
  using navigator_t = caching_navigator<detector_t>;

  // Propagator with pathlimit aborter
  using material_tracer_t = detray::validation::material_tracer<vecmem::vector>;
  using pathlimit_aborter_t = actor::pathlimit_aborter<scalar_t>;
  using parameter_updater_t =
      actor::parameter_updater<algebra_t,
                               actor::pointwise_material_interactor<algebra_t>,
                               material_tracer_t>;
  using actor_chain_t = actor_chain<pathlimit_aborter_t, parameter_updater_t>;
  using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

  // Propagator
  propagator_t prop{cfg};

  // Build actor and propagator states
  typename pathlimit_aborter_t::state pathlimit_aborter_state{
      cfg.stepping.path_limit};
  actor::parameter_updater_state<algebra_t> updater_state{cfg, track};
  typename actor::pointwise_material_interactor<algebra_t>::state
      interactor_state{};
  typename material_tracer_t::state mat_tracer_state{*host_mr};

  auto actor_states = detray::tie(pathlimit_aborter_state, updater_state,
                                  interactor_state, mat_tracer_state);

  typename propagator_t::state propagation{track, det, cfg.context};

  // Run the propagation
  bool success = prop.propagate(propagation, actor_states);

  return std::make_tuple(success,
                         std::move(mat_tracer_state).release_track_material(),
                         std::move(mat_tracer_state).release_material_steps());
}

/// Functor that runs the material recording on the host for a track collection
struct material_recorder {
  /// Test name suffix
  static constexpr std::string_view name{"cpu"};

  /// Run the propagation that records the material
  /// @param host_mr memory resource for host allocations
  /// @param det the detector for which the material should be recorded
  /// @param cfg propagation configuration
  /// @param tracks input tracks for which the material should be collected
  /// @returns the recorded material trace data
  template <concepts::detector detector_t>
  auto operator()(
      vecmem::memory_resource *host_mr, vecmem::memory_resource * /*mr*/,
      const detector_t &det, const propagation::config &cfg,
      const dvector<free_track_parameters<typename detector_t::algebra_type>>
          &tracks,
      const std::vector<std::size_t> & /*capacities*/) const {
    // Prepare output data
    material_trace_data<typename detector_t::algebra_type> mat_traces{};
    mat_traces.tracks.reserve(tracks.size());
    mat_traces.track_materials.reserve(tracks.size());
    mat_traces.material_steps.reserve(tracks.size());

    for (const auto &[i, track] : detray::views::enumerate(tracks)) {
      auto [success, track_mat, mat_steps] =
          detray::validation::host::record_material(host_mr, det, cfg, track);

      if (success) {
        mat_traces.tracks.push_back(track);
        mat_traces.track_materials.push_back(track_mat);
        mat_traces.material_steps.push_back(std::move(mat_steps));
      } else {
        DETRAY_ERROR_HOST("Propagation failed for track "
                          << i << ": Did not collect material trace!");
        DETRAY_DEBUG_HOST("Failed track parameters: " << track);
      }
    }

    return mat_traces;
  }
};

}  // namespace detray::validation::host
