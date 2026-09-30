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
#include "traccc/finding/device/build_tracks.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::build_tracks
__global__ void build_tracks(const bool run_mbf,
                             const measurement_selector::config calib_cfg,
                             const device::build_tracks_payload payload) {
  device::build_tracks(details::global_index1(), run_mbf, calib_cfg, payload);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::build_tracks_kernel(
    unsigned int n_threads, bool run_mbf_smoother,
    const measurement_selector::config& calib_cfg,
    const device::build_tracks_payload& payload) const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = warp_size() * 2;
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;

  // Launch the kernel.
  hipLaunchKernelGGL(kernels::build_tracks, deviceBlocks, deviceThreads, 0,
                     details::get_stream(stream()), run_mbf_smoother, calib_cfg,
                     payload);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

}  // namespace traccc::hip
