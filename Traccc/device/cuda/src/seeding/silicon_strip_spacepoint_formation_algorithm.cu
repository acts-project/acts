// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../sanity/ordered_on.cuh"
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/silicon_strip_spacepoint_formation_algorithm.hpp"

// Project include(s).
#include <thrust/execution_policy.h>
#include <thrust/scan.h>

#include "traccc/seeding/device/count_strip_pairs.hpp"
#include "traccc/seeding/device/find_strip_pairs.hpp"
#include "traccc/seeding/device/form_strip_spacepoints_from_pairs.hpp"
#include "traccc/utils/stream_synchronizing_allocator.hpp"

namespace traccc::cuda {
namespace kernels {

__global__ void count_strip_pairs_kernel(
    edm::measurement_collection::const_view measurements,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    strip_pairing_rule_collection_types::const_view pairing_rules,
    point3 beam_spot, vecmem::data::vector_view<unsigned int> standard_counts,
    vecmem::data::vector_view<unsigned int> overlap_counts) {
  device::count_strip_pairs(details::global_index1(), measurements,
                            surface_infos, pairing_rules, beam_spot,
                            standard_counts, overlap_counts);
}

__global__ void find_strip_pairs_kernel(
    edm::measurement_collection::const_view measurements,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    strip_pairing_rule_collection_types::const_view pairing_rules,
    point3 beam_spot,
    vecmem::data::vector_view<const unsigned int> standard_offsets,
    vecmem::data::vector_view<const unsigned int> overlap_offsets,
    strip_pair_collection_types::view opposite_pairs,
    strip_pair_collection_types::view overlap_pairs) {
  device::find_strip_pairs(details::global_index1(), measurements,
                           surface_infos, pairing_rules, beam_spot,
                           standard_offsets, overlap_offsets, opposite_pairs,
                           overlap_pairs);
}

__global__ void form_strip_spacepoints_from_pairs_kernel(
    edm::measurement_collection::const_view measurements,
    strip_pair_collection_types::const_view pairs,
    strip_measurement_surface_info_collection_types::const_view surface_infos,
    point3 beam_spot, vecmem::data::vector_view<unsigned int> accepted,
    edm::spacepoint_collection::view spacepoints) {
  device::form_strip_spacepoints_from_pairs(details::global_index1(),
                                            measurements, pairs, surface_infos,
                                            beam_spot, accepted, spacepoints);
}

__global__ void gather_strip_spacepoints_kernel(
    edm::spacepoint_collection::const_view candidates,
    vecmem::data::vector_view<const unsigned int> offsets,
    edm::spacepoint_collection::view spacepoints) {
  device::gather_strip_spacepoints(details::global_index1(), candidates,
                                   offsets, spacepoints);
}

struct measurement_surface_order {
  template <typename T1, typename T2>
  TRACCC_HOST_DEVICE bool operator()(const edm::measurement<T1>& a,
                                     const edm::measurement<T2>& b) const {
    return a.surface_link().index() <= b.surface_link().index();
  }
};

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
      payload.beam_spot, payload.standard_counts, payload.overlap_counts);
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
      payload.beam_spot, payload.standard_offsets, payload.overlap_offsets,
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
      payload.beam_spot, payload.accepted, payload.spacepoints);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
  stream().synchronize();
}

void silicon_strip_spacepoint_formation_algorithm::scan_offsets(
    vecmem::data::vector_view<unsigned int> offsets) const {
  assert(offsets.size_ptr() == nullptr);
  thrust::inclusive_scan(
      thrust::cuda::par_nosync(
          stream_synchronizing_allocator(mr().main, stream()))
          .on(details::get_stream(stream())),
      offsets.ptr(), offsets.ptr() + offsets.capacity(), offsets.ptr());
}

bool silicon_strip_spacepoint_formation_algorithm::input_is_sorted(
    const edm::measurement_collection::const_view& measurements) const {
  return is_ordered_on<edm::measurement_collection::const_device>(
      kernels::measurement_surface_order{}, mr().main, copy(), stream(),
      measurements);
}

void silicon_strip_spacepoint_formation_algorithm::gather_spacepoints_kernel(
    const gather_spacepoints_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (payload.n_pairs + n_threads - 1u) / n_threads;
  kernels::gather_strip_spacepoints_kernel<<<n_blocks, n_threads, 0,
                                             details::get_stream(stream())>>>(
      payload.candidates, payload.offsets, payload.spacepoints);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
  // The base class destroys candidate and flag buffers after this call.
  stream().synchronize();
}

}  // namespace traccc::cuda
