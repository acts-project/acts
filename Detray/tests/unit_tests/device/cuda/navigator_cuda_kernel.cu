// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/definitions/detail/cuda_definitions.hpp"

// Detray test include(s)
#include "navigator_cuda_kernel.hpp"

namespace detray {

__global__ void navigator_test_kernel(
    typename detector_host_t::view_type det_data,
    const propagation::config prop_cfg,
    vecmem::data::vector_view<free_track_parameters<test_algebra>> tracks_data,
    vecmem::data::jagged_vector_view<dindex> volume_records_data,
    vecmem::data::jagged_vector_view<point3> position_records_data) {
  int gid = threadIdx.x + blockIdx.x * blockDim.x;

  detector_device_t det(det_data);
  vecmem::device_vector<free_track_parameters<test_algebra>> tracks(
      tracks_data);
  vecmem::jagged_device_vector<dindex> volume_records(volume_records_data);
  vecmem::jagged_device_vector<point3> position_records(position_records_data);

  if (gid >= tracks.size()) {
    return;
  }

  navigator_device_t nav;

  auto& traj = tracks.at(gid);
  stepper_t stepper;

  const navigation::config& nav_cfg = prop_cfg.navigation;
  const stepping::config& step_cfg = prop_cfg.stepping;

  prop_state<navigator_device_t::state> propagation{
      stepper_t::state{traj}, navigator_device_t::state(det)};

  navigator_device_t::state& navigation = propagation.navigation();
  stepper_t::state& stepping = propagation.stepping();
  const auto& ctx = propagation.context();

  // Set initial volume
  navigation.set_volume(0u);

  // Start propagation and record volume IDs
  nav.init(stepping(), navigation, nav_cfg, ctx);
  bool heartbeat = navigation.is_alive();
  bool do_reset{true};

  while (heartbeat) {
    heartbeat &= stepper.step(navigation(), stepping, step_cfg, do_reset);

    navigation.set_high_trust();

    do_reset = nav.update(stepping(), navigation, nav_cfg, prop_cfg.context);
    do_reset |= navigation.is_on_surface();
    heartbeat &= navigation.is_alive();

    // Record volume
    volume_records[gid].push_back(navigation.volume());
    position_records[gid].push_back(stepping().pos());
  }
}

void navigator_test(
    typename detector_host_t::view_type det_data, propagation::config& prop_cfg,
    vecmem::data::vector_view<free_track_parameters<test_algebra>>& tracks_data,
    vecmem::data::jagged_vector_view<dindex>& volume_records_data,
    vecmem::data::jagged_vector_view<point3>& position_records_data) {
  constexpr int thread_dim = 2 * WARP_SIZE;
  constexpr int block_dim = theta_steps * phi_steps / thread_dim + 1;

  // run the test kernel
  navigator_test_kernel<<<block_dim, thread_dim>>>(
      det_data, prop_cfg, tracks_data, volume_records_data,
      position_records_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
