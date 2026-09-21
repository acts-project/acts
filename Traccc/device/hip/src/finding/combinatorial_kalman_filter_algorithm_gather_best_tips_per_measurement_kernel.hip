/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/barrier.hpp"
#include "../utils/global_index.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/finding/device/gather_best_tips_per_measurement.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running
/// @c traccc::device::gather_best_tips_per_measurement
__global__ void gather_best_tips_per_measurement(
    const device::gather_best_tips_per_measurement_payload<default_algebra>
        payload) {
  device::gather_best_tips_per_measurement(details::global_index1(), barrier{},
                                           payload);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::
    gather_best_tips_per_measurement_kernel(
        unsigned int n_threads,
        const device::gather_best_tips_per_measurement_payload<default_algebra>&
            payload) const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = warp_size();
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;

  // Launch the kernel.
  hipLaunchKernelGGL(kernels::gather_best_tips_per_measurement, deviceBlocks,
                     deviceThreads, 0, details::get_stream(stream()), payload);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

}  // namespace traccc::hip
