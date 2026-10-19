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
  std::cout << "We are inside the hough_transform_seeding_algorithm::operator()" << std::endl;
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

  std::cout << "Hough Transform will run on a total of " << n_spacepoints << " input space points" << std::endl;
  // If there are no spacepoints, return right away.
  if (n_spacepoints == 0) {
    std::cout << "well, not really ..." << std::endl;
    return {};
  }

  // Set up the output buffer.
  edm::houghplane_collection<N_COLS, N_ROWS>::const_view::size_type n_slices = 25u;
  edm::houghplane_collection<N_COLS, N_ROWS>::buffer result(n_slices, mr().main);
  copy().setup(result)->ignore();

  // Fill plane edges
  

  // Fill the planes with space points
  //
  //
  //
  // Launch the seed parameter estimation kernel.
  // estimate_seed_params_kernel({n_seeds, m_data->m_config, bfield, measurements,
  //                             spacepoints, seeds, result});

  // Return the result.
  return result;
}

}  // namespace traccc::device
