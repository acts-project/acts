// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Library include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/hough_transform_seeding_algorithm.hpp"

namespace traccc::cuda {
namespace kernels {
   __global__ 
   void define_plane_edges(const traccc::hough_seeding_config config,
		   edm::houghplane_collection<device::N_COLS, device::N_ROWS>::view houghplanes) {
      float steps_x = (config.max_range_x - config.min_range_x) / device::N_COLS;
      float steps_y = (config.max_range_y - config.min_range_y) / device::N_ROWS;
      device::define_plane_edges(details::global_index1(), houghplanes,
		     steps_x, config.min_range_x,
		     steps_y, config.min_range_y); 
   }
      
   __global__
   void fill_hough_plane(const traccc::hough_seeding_config config,
		   edm::spacepoint_collection::const_view::size_type n_spacepoints,
		   const edm::spacepoint_collection::const_view& spacepoints,
		   edm::houghplane_collection<device::N_COLS, device::N_ROWS>::view& houghplanes) {
      device::fill_hough_plane(details::global_index1(), n_spacepoints, spacepoints, 
		      config.interaction_center, config.interaction_window,
		      houghplanes);
   }
} // namespace kernels

hough_transform_seeding_algorithm::hough_transform_seeding_algorithm(
    const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    const hough_seeding_config& seeding_config,
    std::unique_ptr<const Logger> logger, await_function_type await_func)
    : device::hough_transform_seeding_algorithm(mr, copy, seeding_config,
                                        std::move(logger)),
      cuda::algorithm_base{str, std::move(await_func)} {}

void hough_transform_seeding_algorithm::define_plane_edges_kernel(const define_plane_edges_kernel_payload& payload) const {
   unsigned int maxValues = std::max(device::N_COLS, device::N_ROWS);
   const unsigned int n_threads = warp_size() * 4;
   const unsigned int n_blocks = (maxValues  + n_threads - 1) / n_threads;
   kernels::define_plane_edges<<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config, payload.houghplanes);
}

void hough_transform_seeding_algorithm::fill_hough_plane_kernel(const fill_hough_planes_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 4;
  const unsigned int n_blocks = (payload.n_spacepoints + n_threads - 1) / n_threads;
  kernels::fill_hough_plane<<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config, payload.n_spacepoints, payload.spacepoints, payload.houghplanes);
}

void hough_transform_seeding_algorithm::synchronize() const {
  stream().synchronize();
}

}  // namespace traccc::cuda
