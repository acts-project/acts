// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/silicon_strip_spacepoint_formation_algorithm.hpp"

// Project include(s).
#include "traccc/seeding/device/count_strip_pairs.hpp"
#include "traccc/seeding/device/find_strip_pairs.hpp"
#include "traccc/seeding/device/form_strip_spacepoints_from_pairs.hpp"

namespace traccc::cuda {
namespace kernels {

__global__ void __launch_bounds__(1024, 1) count_strip_pairs_kernel(
    edm::measurement_collection::const_view measurements,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    strip_pairing_rule_collection_types::const_view pairing_rules,
    point3 beam_spot, unsigned int& n_opposite_pairs,
    unsigned int& n_overlap_pairs) {
  device::count_strip_pairs(details::global_index1(), measurements,
                            surface_infos, pairing_rules, beam_spot,
                            n_opposite_pairs, n_overlap_pairs);
}

__global__ void __launch_bounds__(1024, 1) find_strip_pairs_kernel(
    edm::measurement_collection::const_view measurements,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    strip_pairing_rule_collection_types::const_view pairing_rules,
    point3 beam_spot, unsigned int& opposite_position,
    unsigned int& overlap_position,
    strip_pair_collection_types::view opposite_pairs,
    strip_pair_collection_types::view overlap_pairs) {
  device::find_strip_pairs(details::global_index1(), measurements,
                           surface_infos, pairing_rules, beam_spot,
                           opposite_position, overlap_position, opposite_pairs,
                           overlap_pairs);
}

__global__ void
__launch_bounds__(1024, 1) form_strip_spacepoints_from_pairs_kernel(
    edm::measurement_collection::const_view measurements,
    strip_pair_collection_types::const_view pairs,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    point3 beam_spot, edm::spacepoint_collection::view spacepoints) {
  device::form_strip_spacepoints_from_pairs(details::global_index1(),
                                            measurements, pairs, surface_infos,
                                            beam_spot, spacepoints);
}

}  // namespace kernels

silicon_strip_spacepoint_formation_algorithm::
    silicon_strip_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, const vecmem::copy& copy,
        const stream_wrapper& str, std::unique_ptr<const Logger> logger)
    : device::silicon_strip_spacepoint_formation_algorithm(mr, copy,
                                                           std::move(logger)),
      cuda::algorithm_base(str) {}

void silicon_strip_spacepoint_formation_algorithm::count_strip_pairs_kernel(
    const count_strip_pairs_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks =
      (payload.n_measurements + n_threads - 1) / n_threads;
  kernels::count_strip_pairs_kernel<<<n_blocks, n_threads, 0,
                                      details::get_stream(stream())>>>(
      payload.measurements, payload.surface_infos, payload.pairing_rules,
      payload.beam_spot, payload.n_opposite_pairs, payload.n_overlap_pairs);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void silicon_strip_spacepoint_formation_algorithm::find_strip_pairs_kernel(
    const find_strip_pairs_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks =
      (payload.n_measurements + n_threads - 1) / n_threads;
  kernels::find_strip_pairs_kernel<<<n_blocks, n_threads, 0,
                                     details::get_stream(stream())>>>(
      payload.measurements, payload.surface_infos, payload.pairing_rules,
      payload.beam_spot, payload.opposite_position, payload.overlap_position,
      payload.opposite_pairs, payload.overlap_pairs);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void silicon_strip_spacepoint_formation_algorithm::form_spacepoints_kernel(
    const form_spacepoints_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (payload.n_pairs + n_threads - 1) / n_threads;
  kernels::form_strip_spacepoints_from_pairs_kernel<<<
      n_blocks, n_threads, 0, details::get_stream(stream())>>>(
      payload.measurements, payload.pairs, payload.surface_infos,
      payload.beam_spot, payload.spacepoints);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
  stream().synchronize();
}

}  // namespace traccc::cuda
