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
#include "traccc/fitting/device/fill_fitting_sort_keys.hpp"

// (Roc)Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/sort.h>

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <cassert>
#include <memory_resource>

namespace traccc::hip {
namespace kernels {

__global__ void fill_fitting_sort_keys(
    edm::track_collection<default_algebra>::const_view track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view) {
  device::fill_fitting_sort_keys(details::global_index1(),
                                 track_candidates_view, keys_view, ids_view);
}

}  // namespace kernels

void kalman_fitting_algorithm::prepare_track_fit_order(
    const edm::track_collection<default_algebra>::const_view& tracks,
    vecmem::data::vector_view<device::sort_key>& track_sort_keys,
    vecmem::data::vector_view<unsigned int>& track_indices) const {
  // Get the number of tracks.
  const unsigned int n_tracks = tracks.capacity();
  assert(n_tracks == copy().get_size(tracks));
  assert(n_tracks == track_indices.capacity());
  assert(track_indices.size_ptr() == nullptr);

  // Launch parameters for the kernel.
  const unsigned int nThreads = warp_size() * 4;
  const unsigned int nBlocks = (n_tracks + nThreads - 1) / nThreads;

  // Get a convenience variable for the stream that we'll be using.
  hipStream_t str = details::get_stream(stream());

  // Fill the keys and indices buffers.
  hipLaunchKernelGGL(kernels::fill_fitting_sort_keys, nBlocks, nThreads, 0, str,
                     tracks, track_sort_keys, track_indices);

// Set up the Thrust execution policy. Unfortunately this has to be done in
// a backend specify way. :-(
#if defined(__HIP_PLATFORM_AMD__)
  auto policy =
      thrust::hip::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(str);
#elif defined(__HIP_PLATFORM_NVIDIA__)
  auto policy =
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(str);
#else
#error "Unrecognised HIP platform."
#endif

  // Sort the key to get the sorted parameter ids
  vecmem::device_vector<device::sort_key> keys_device(track_sort_keys);
  vecmem::device_vector<unsigned int> track_indices_device(track_indices);
  thrust::sort_by_key(policy, keys_device.begin(), keys_device.end(),
                      track_indices_device.begin());
}

}  // namespace traccc::hip
