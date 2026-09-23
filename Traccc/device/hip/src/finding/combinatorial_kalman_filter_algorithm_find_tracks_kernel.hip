/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/barrier.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/finding/device/find_tracks.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <utility>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::find_tracks
///
/// The kernel is not constrained on @c detray::concepts::detector, because
/// hip-clang checks the constraint in the device context of the kernel,
/// where the host-only @c name_map constructor that the concept uses is not
/// callable. The device function that it calls carries the constraint.
///
template <detray::concepts::detector detector_t>
__global__ void find_tracks(const finding_config cfg,
                            const typename detector_t::const_view_type det,
                            const device::find_tracks_payload payload) {
  __shared__ unsigned int shared_num_out_params;
  __shared__ unsigned int shared_candidates_size;
  extern __shared__ unsigned long long int s[];
  unsigned long long int* shared_insertion_mutex = s;
  std::pair<unsigned int, unsigned int>* shared_candidates =
      reinterpret_cast<std::pair<unsigned int, unsigned int>*>(
          &shared_insertion_mutex[blockDim.x]);

  hip::barrier barrier;
  details::thread_id1 thread_id;

  device::find_tracks<detector_t>(
      thread_id, barrier, cfg, det, payload,
      device::find_tracks_shared_payload{
          .shared_num_out_params = shared_num_out_params,
          .shared_insertion_mutex = shared_insertion_mutex,
          .shared_candidates = shared_candidates,
          .shared_candidates_size = shared_candidates_size});
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::find_tracks_kernel(
    unsigned int n_threads, const finding_config& config,
    const detector_buffer& detector,
    const device::find_tracks_payload& payload) const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = warp_size() * 2;
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;
  const std::size_t deviceSharedMem =
      deviceThreads * sizeof(unsigned long long int) +
      2 * deviceThreads * sizeof(std::pair<unsigned int, unsigned int>);

  // Launch the kernel for the appropriate detector type.
  detector_buffer_visitor<detector_type_list>(
      detector, [&]<detray::concepts::detector detector_t>(
                    const detray::detector_view_t<detector_t>& det) {
        using detector_device_t = detray::detector_device_t<detector_t>;

        hipLaunchKernelGGL(kernels::find_tracks<detector_device_t>,
                           deviceBlocks, deviceThreads, deviceSharedMem,
                           details::get_stream(stream()), config, det, payload);
        TRACCC_HIP_ERROR_CHECK(hipGetLastError());
      });
}

}  // namespace traccc::hip
