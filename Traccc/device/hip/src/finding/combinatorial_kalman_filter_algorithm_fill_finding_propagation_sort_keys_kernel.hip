/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/global_index.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/finding/device/fill_finding_propagation_sort_keys.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::fill_finding_propagation_sort_keys
__global__ void fill_finding_propagation_sort_keys(
    const device::fill_finding_propagation_sort_keys_payload payload) {
  device::fill_finding_propagation_sort_keys(details::global_index1(), payload);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::
    fill_finding_propagation_sort_keys_kernel(
        unsigned int n_threads,
        const device::fill_finding_propagation_sort_keys_payload& payload)
        const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = 64;
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;

  // Launch the kernel.
  hipLaunchKernelGGL(kernels::fill_finding_propagation_sort_keys, deviceBlocks,
                     deviceThreads, 0, details::get_stream(stream()), payload);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

}  // namespace traccc::hip
