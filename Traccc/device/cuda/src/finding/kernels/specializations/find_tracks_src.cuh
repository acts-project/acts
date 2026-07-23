/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/finding/device/find_tracks.hpp"

// Local include(s).
#include "../../../utils/barrier.hpp"
#include "../../../utils/hints.hpp"
#include "../../../utils/thread_id.hpp"

// System include(s).
#include <utility>

namespace traccc::cuda {
namespace kernels {

template <typename detector_t>
__global__ void find_tracks(
    const __grid_constant__ finding_config cfg,
    const __grid_constant__ typename detector_t::const_view_type det,
    const __grid_constant__ device::find_tracks_payload payload) {
  __shared__ unsigned int shared_num_out_params;
  __shared__ unsigned int shared_candidates_size;
  extern __shared__ unsigned long long int s[];
  unsigned long long int* shared_insertion_mutex = s;
  std::pair<unsigned int, unsigned int>* shared_candidates =
      reinterpret_cast<std::pair<unsigned int, unsigned int>*>(
          &shared_insertion_mutex[blockDim.x]);

  cuda::barrier barrier;
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

template <typename detector_t>
void find_tracks(const dim3& grid_size, const dim3& block_size,
                 std::size_t shared_mem_size, const cudaStream_t& stream,
                 const finding_config& cfg,
                 const typename detector_t::const_view_type& det,
                 const device::find_tracks_payload& payload) {
  kernels::find_tracks<detector_t>
      <<<grid_size, block_size, shared_mem_size, stream>>>(cfg, det, payload);
}
}  // namespace traccc::cuda
