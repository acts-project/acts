// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/definitions/detail/cuda_definitions.hpp"
#include "propagation.hpp"

namespace detray::tutorial {

// Propagation configurations
inline constexpr scalar path_limit{2.f * detray::unit<scalar>::m};

/// Kernel that runs the entire propagation loop
__global__ void propagation_kernel(
    typename detray::tutorial::detector_host_t::view_type det_data,
    typename detray::tutorial::device_field_t::view_t field_data,
    const vecmem::data::vector_view<detray::tutorial::track_t> tracks_data) {
  int gid = threadIdx.x + blockIdx.x * blockDim.x;

  // Setup device-side track collection
  vecmem::device_vector<detray::tutorial::track_t> tracks(tracks_data);

  if (gid >= tracks.size()) {
    return;
  }

  // Setup of the device-side detector
  detray::tutorial::detector_device_t det(det_data);

  // Create propagator from a stepper and a navigator
  propagation::config cfg{};
  cfg.navigation.search_window = {3u, 3u};
  detray::tutorial::propagator_t p{cfg};

  // Create actor states
  detray::actor::pathlimit_aborter<scalar>::state aborter_state{path_limit};
  detray::actor::parameter_updater_state<algebra_t> updater_state{cfg};
  detray::actor::pointwise_material_interactor<algebra_t>::state
      interactor_state{};

  auto actor_states =
      detray::tie(aborter_state, updater_state, interactor_state);

  // Create the propagator state for the track
  detray::tutorial::propagator_t::state state(tracks[gid], field_data, det);

  // Run propagation
  p.propagate(state, actor_states);
}

void propagation(typename detray::tutorial::detector_host_t::view_type det_data,
                 typename detray::tutorial::device_field_t::view_t field_data,
                 const vecmem::data::vector_view<track_t> tracks_data) {
  int thread_dim = 2 * WARP_SIZE;
  int block_dim = tracks_data.size() / thread_dim + 1;

  // run the tutorial kernel
  propagation_kernel<<<block_dim, thread_dim>>>(det_data, field_data,
                                                tracks_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray::tutorial
