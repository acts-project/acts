// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "navigator_cuda_kernel.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"

// vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

TEST(navigator_cuda, navigator) {
  // Helper object for performing memory copies.
  vecmem::cuda::copy copy;

  // VecMem memory resource(s)
  vecmem::cuda::managed_memory_resource mng_mr;
  vecmem::cuda::device_memory_resource dev_mr;

  // Create detector
  auto [det, names] = build_toy_detector<test_algebra>(mng_mr);

  propagation::config prop_cfg{};
  navigation::config& nav_cfg = prop_cfg.navigation;
  stepping::config& step_cfg = prop_cfg.stepping;

  // Create navigator
  navigator_host_t nav;
  nav_cfg.search_window = {3u, 3u};

  // Create the vector of initial track parameters
  vecmem::vector<free_track_parameters<test_algebra>> tracks_host(&mng_mr);
  vecmem::vector<free_track_parameters<test_algebra>> tracks_device(&mng_mr);

  // Magnitude of total momentum of tracks
  const scalar p_mag{10.f * unit<scalar>::GeV};

  // Iterate through uniformly distributed momentum directions
  for (auto track :
       uniform_track_generator<free_track_parameters<test_algebra>>(
           phi_steps, theta_steps, p_mag)) {
    tracks_host.push_back(track);
    tracks_device.push_back(track);
  }

  // Host volume index and track position records
  vecmem::jagged_vector<dindex> volume_records_host(theta_steps * phi_steps,
                                                    &mng_mr);
  vecmem::jagged_vector<point3> position_records_host(theta_steps * phi_steps,
                                                      &mng_mr);

  for (unsigned int i = 0u; i < theta_steps * phi_steps; i++) {
    auto& track = tracks_host[i];
    stepper_t stepper;

    // Propagator is built from the stepper and navigator
    prop_state<navigator_host_t::state> propagation{
        stepper_t::state{track}, navigator_host_t::state(det)};

    navigator_host_t::state& navigation = propagation.navigation();
    stepper_t::state& stepping = propagation.stepping();
    const auto& ctx = propagation.context();

    // Start propagation and record volume IDs
    nav.init(stepping(), navigation, nav_cfg, ctx);
    bool heartbeat = navigation.is_alive();
    bool do_reset{true};

    while (heartbeat) {
      heartbeat =
          heartbeat && stepper.step(navigation(), stepping, step_cfg, do_reset);

      navigation.set_high_trust();

      do_reset = nav.update(stepping(), navigation, nav_cfg, prop_cfg.context);
      do_reset = do_reset || navigation.is_on_surface();
      heartbeat = heartbeat && navigation.is_alive();

      // Record volume
      volume_records_host[i].push_back(navigation.volume());
      position_records_host[i].push_back(stepping().pos());
    }
  }

  // Device volume index and track position records
  vecmem::jagged_vector<dindex> volume_records_device(&mng_mr);
  vecmem::jagged_vector<point3> position_records_device(&mng_mr);

  // Create size and capacity vectors for volume record buffer
  std::vector<std::size_t> capacities;

  for (unsigned int i = 0u; i < theta_steps * phi_steps; i++) {
    capacities.push_back(volume_records_host[i].size());
  }

  vecmem::data::jagged_vector_buffer<dindex> volume_records_buffer(
      capacities, dev_mr, &mng_mr, vecmem::data::buffer_type::resizable);
  copy.setup(volume_records_buffer)->wait();

  vecmem::data::jagged_vector_buffer<point3> position_records_buffer(
      capacities, dev_mr, &mng_mr, vecmem::data::buffer_type::resizable);
  copy.setup(position_records_buffer)->wait();

  // Get detector data
  auto det_data = detray::get_data(det);

  // Get tracks data
  auto tracks_data = vecmem::get_data(tracks_device);

  // Run navigator test
  navigator_test(det_data, prop_cfg, tracks_data, volume_records_buffer,
                 position_records_buffer);

  // Copy volume record buffer into volume & position records device
  copy(volume_records_buffer, volume_records_device)->wait();
  copy(position_records_buffer, position_records_device)->wait();

  for (unsigned int i = 0u; i < volume_records_host.size(); i++) {
    EXPECT_EQ(volume_records_host[i].size(), volume_records_device[i].size());

    for (unsigned int j = 0u; j < volume_records_host[i].size(); j++) {
      EXPECT_EQ(volume_records_host[i][j], volume_records_device[i][j]);

      auto& pos_host = position_records_host[i][j];
      auto& pos_device = position_records_device[i][j];

      EXPECT_NEAR(pos_host[0], pos_device[0], pos_diff_tolerance);
      EXPECT_NEAR(pos_host[1], pos_device[1], pos_diff_tolerance);
      EXPECT_NEAR(pos_host[2], pos_device[2], pos_diff_tolerance);
    }
  }
}
