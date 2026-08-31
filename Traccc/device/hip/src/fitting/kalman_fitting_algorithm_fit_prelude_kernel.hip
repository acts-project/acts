/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/global_index.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/fitting/kalman_fitting_algorithm.hpp"

// Project include(s).
#include "traccc/fitting/device/fit_prelude.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <cassert>

namespace traccc::hip {
namespace kernels {

__global__ void fit_prelude(const device::fit_prelude_payload payload) {
  device::fit_prelude(details::global_index1(), payload);
}

}  // namespace kernels

void kalman_fitting_algorithm::fit_prelude_kernel(
    const device::fit_prelude_payload& payload) const {
  // Get the number of tracks.
  const unsigned int n_tracks = payload.input_tracks.tracks.capacity();
  assert(n_tracks == copy().get_size(payload.input_tracks.tracks));
  assert(n_tracks == payload.track_indices.capacity());
  assert(payload.track_indices.size_ptr() == nullptr);
  assert(n_tracks == copy().get_size(payload.output_tracks.tracks));

  // Launch parameters for the kernel.
  const unsigned int nThreads = warp_size() * 4;
  const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

  // Run the fitting, using the sorted parameter IDs.
  hipLaunchKernelGGL(kernels::fit_prelude, nBlocks, nThreads, 0,
                     details::get_stream(stream()), payload);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

}  // namespace traccc::hip
