// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Library include(s).
#include "../utils/hip_error_handling.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/clusterization/measurement_sorting_algorithm.hpp"

// Project include(s).
#include "traccc/clusterization/device/measurement_sorting.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>

// Project include(s).
#include "traccc/utils/stream_synchronizing_allocator.hpp"

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/sort.h>

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// Kernel wrapping @c traccc::device::fill_measurement_sort_keys
__global__ void fill_measurement_sort_keys(
    const edm::measurement_collection::const_view measurements_view,
    vecmem::data::vector_view<device::measurement_sort_key_t> keys_view,
    vecmem::data::vector_view<unsigned int> indices_view) {
  device::fill_measurement_sort_keys(details::thread_id1{}.getGlobalThreadId(),
                                     measurements_view, keys_view,
                                     indices_view);
}

/// Kernel wrapping @c traccc::device::fill_sorted_measurements
__global__ void fill_sorted_measurements(
    const edm::measurement_collection::const_view input_view,
    edm::measurement_collection::view output_view,
    const vecmem::data::vector_view<const unsigned int> sorted_indices_view) {
  device::fill_sorted_measurements(details::thread_id1{}.getGlobalThreadId(),
                                   input_view, output_view,
                                   sorted_indices_view);
}

}  // namespace kernels

measurement_sorting_algorithm::measurement_sorting_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    const stream_wrapper& str, std::unique_ptr<const Logger> logger)
    : device::measurement_sorting_algorithm(mr, copy, std::move(logger)),
      hip::algorithm_base(str) {}

void measurement_sorting_algorithm::sorting_kernel(
    const measurement_sorting_kernel_payload& payload) const {
  const unsigned int n_measurements = payload.n_measurements;
  hipStream_t str = details::get_stream(stream());
  // Set up the Thrust execution policy. Unfortunately this has to be done in
  // a backend specify way. :-(
#if defined(__HIP_PLATFORM_AMD__)
  auto policy = thrust::hip::par_nosync(
                    stream_synchronizing_allocator(mr().main, stream()))
                    .on(str);
#elif defined(__HIP_PLATFORM_NVIDIA__)
  auto policy = thrust::cuda::par_nosync(
                    stream_synchronizing_allocator(mr().main, stream()))
                    .on(str);
#else
#error "Unrecognised HIP platform."
#endif

  const unsigned int num_threads = warp_size() * 8;
  const unsigned int num_blocks =
      (n_measurements + num_threads - 1) / num_threads;

  // Sort the indices by the sorting keys, with a radix sort.
  hipLaunchKernelGGL(kernels::fill_measurement_sort_keys, num_blocks,
                     num_threads, 0, str, payload.measurements, payload.keys,
                     payload.indices);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
  thrust::sort_by_key(policy, payload.keys.ptr(),
                      payload.keys.ptr() + n_measurements,
                      payload.indices.ptr());

  // Fill the output with the sorted measurements.
  hipLaunchKernelGGL(kernels::fill_sorted_measurements, num_blocks, num_threads,
                     0, str, payload.measurements, payload.output,
                     payload.indices);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

void measurement_sorting_algorithm::synchronize() const {
  stream().synchronize();
}

}  // namespace traccc::hip
