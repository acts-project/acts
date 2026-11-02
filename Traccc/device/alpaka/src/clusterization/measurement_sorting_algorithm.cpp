// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/alpaka/clusterization/measurement_sorting_algorithm.hpp"

#include "../utils/get_queue.hpp"
#include "../utils/parallel_algorithms.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/clusterization/device/measurement_sorting.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>

// System include(s).
#include <memory_resource>

namespace traccc::alpaka {
namespace kernels {

/// Kernel wrapping @c traccc::device::fill_measurement_sort_keys
struct fill_measurement_sort_keys {
  /// @param[in]  acc               Alpaka accelerator object
  /// @param[in]  measurements_view View of the unsorted measurements
  /// @param[out] keys_view         View of the sorting keys
  /// @param[out] indices_view      View of the measurement indices
  ///
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const edm::measurement_collection::const_view measurements_view,
      vecmem::data::vector_view<device::measurement_sort_key_t> keys_view,
      vecmem::data::vector_view<unsigned int> indices_view) const {
    device::fill_measurement_sort_keys(
        details::thread_id1{acc}.getGlobalThreadId(), measurements_view,
        keys_view, indices_view);
  }
};  // struct fill_measurement_sort_keys

/// Kernel filling the output buffer with sorted measurements.
struct fill_sorted_measurements {
  /// @param[in] acc Alpaka accelerator object
  /// @param[in] input_view View of the input measurements
  /// @param[out] output_view View of the output measurements
  /// @param[in] sorted_indices_view View of the sorted measurement indices
  ///
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc, const edm::measurement_collection::const_view input_view,
      edm::measurement_collection::view output_view,
      const vecmem::data::vector_view<const unsigned int> sorted_indices_view)
      const {
    device::fill_sorted_measurements(
        details::thread_id1{acc}.getGlobalThreadId(), input_view, output_view,
        sorted_indices_view);
  }
};  // struct fill_sorted_measurements

}  // namespace kernels

measurement_sorting_algorithm::measurement_sorting_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    alpaka::queue& q, std::unique_ptr<const Logger> logger,
    await_function_type await_func)
    : device::measurement_sorting_algorithm(mr, copy, std::move(logger)),
      alpaka::algorithm_base(q, std::move(await_func)) {}

void measurement_sorting_algorithm::sorting_kernel(
    const measurement_sorting_kernel_payload& payload) const {
  const unsigned int n_measurements = payload.n_measurements;

  const unsigned int num_threads = warp_size() * 8;
  const unsigned int num_blocks =
      (n_measurements + num_threads - 1) / num_threads;
  auto workDiv = makeWorkDiv<Acc>(num_blocks, num_threads);

  // Sort the indices by the sorting keys, with a radix sort.
  ::alpaka::exec<Acc>(details::get_queue(queue()), workDiv,
                      kernels::fill_measurement_sort_keys{},
                      payload.measurements, vecmem::get_data(payload.keys),
                      vecmem::get_data(payload.indices));
  details::sort_by_key(queue(), mr(), payload.keys.ptr(),
                       payload.keys.ptr() + n_measurements,
                       payload.indices.ptr());

  // Fill the output with the sorted measurements.
  ::alpaka::exec<Acc>(details::get_queue(queue()), workDiv,
                      kernels::fill_sorted_measurements{}, payload.measurements,
                      vecmem::get_data(payload.output),
                      vecmem::get_data(payload.indices));
}

void measurement_sorting_algorithm::synchronize() const {
  queue().synchronize();
}

}  // namespace traccc::alpaka
