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

// Project include(s).
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/seeding/device/count_doublets.hpp"
#include "traccc/seeding/device/count_grid_capacities.hpp"
#include "traccc/seeding/device/count_triplets.hpp"
#include "traccc/seeding/device/find_doublets.hpp"
#include "traccc/seeding/device/find_triplets.hpp"
#include "traccc/seeding/device/populate_grid.hpp"
#include "traccc/seeding/device/reduce_triplet_counts.hpp"
#include "traccc/seeding/device/select_seeds.hpp"
#include "traccc/seeding/device/update_triplet_weights.hpp"

namespace traccc::cuda {
namespace kernels {
   __global__ 
   void define_plane_edges(edm::houghplane_collection<device::N_COLS, device::N_ROWS>::view houghplanes) {
      // we need a config actually
      float minQoptRange = -1.1f;
      float maxQoptRange = 1.1f;
      float minPhi = -std::numbers::pi_v<float>;
      float maxPhi = std::numbers::pi_v<float>;
      float stepQopt = (maxQoptRange - minQoptRange) / device::N_COLS;
      float stepPhi = (maxPhi - minPhi) / device::N_ROWS;
      device::define_plane_edges(details::global_index1(), houghplanes, 
		      stepQopt, minQoptRange,
		      stepPhi, minPhi);
   }
} // namespace kernels

hough_transform_seeding_algorithm::hough_transform_seeding_algorithm(
    const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    std::unique_ptr<const Logger> logger, await_function_type await_func)
    : device::hough_transform_seeding_algorithm(mr, copy,
                                        std::move(logger)),
      cuda::algorithm_base{str, std::move(await_func)} {}

void hough_transform_seeding_algorithm::define_plane_edges_kernel(const define_plane_edges_kernel_payload& payload) const {
   unsigned int maxValues = std::max(device::N_COLS, device::N_ROWS);
   const unsigned int n_threads = warp_size() * 4;
   const unsigned int n_blocks = (maxValues  + n_threads - 1) / n_threads;
   kernels::define_plane_edges<<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.houghplanes);
}

void hough_transform_seeding_algorithm::fill_hough_planes_kernel(const fill_hough_planes_kernel_payload& payload) const {}

void hough_transform_seeding_algorithm::synchronize() const {
  stream().synchronize();
}

}  // namespace traccc::cuda
