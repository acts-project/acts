// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/seeding/device/hough_transform_seeding_algorithm.hpp"

namespace traccc::device {


hough_transform_seeding_algorithm::hough_transform_seeding_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      algorithm_base{mr, copy} {}

auto hough_transform_seeding_algorithm::operator()(
    const edm::spacepoint_collection::const_view& spacepoints) const -> output_type {
      	// Get the number of space points. In an asynchronous way if possible.
  edm::spacepoint_collection::const_view::size_type n_spacepoints = 0u;
  if (mr().host) {
    vecmem::async_size size = copy().get_size(spacepoints, *(mr().host));
    // Block or suspend execution until the size is available.
    await(size);
    n_spacepoints = size.unsafe_get();
  } else {
    n_spacepoints = copy().get_size(spacepoints);
  }

  // If there are no spacepoints, return right away.
  if (n_spacepoints == 0) {
    return {};
  }

  // Set up the output buffer.
  edm::houghplane_collection<N_COLS, N_ROWS>::const_view::size_type n_slices = 25u;
  edm::houghplane_collection<N_COLS, N_ROWS>::buffer result(n_slices, mr().main);
  copy().setup(result)->ignore();

  // Fill plane edges
  // we will move to a config at a certain point and template the algorithm
  define_plane_edges_kernel({device::N_COLS, device::N_ROWS, result});

  /*
  edm::houghplane_collection<N_COLS, N_ROWS>::const_view host_view;
  edm::houghplane_collection<N_COLS, N_ROWS>::buffer host_buffer;
  if (mr().host) {
    host_buffer = copy().to(result, *(mr().host));
    host_view = host_buffer;
  } else {
    host_view = result;
  }

  edm::houghplane_collection<N_COLS, N_ROWS>::const_device host_device{host_view};
  for (unsigned int i = 0; i < host_device.size(); ++i) {
    const edm::houghplane plane = host_device.at(i);
    const auto& edges_x = plane.edges_x();
    const auto& edges_y = plane.edges_y();
  }
  */

  // Fill the planes with space points
  //
  //

  synchronize();

  // Return the result.
  return result;
}

}  // namespace traccc::device
