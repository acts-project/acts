/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/global_index.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/finding/device/condense_tracks.hpp"

// (Roc)Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/scan.h>

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <cassert>
#include <memory_resource>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::condense_tracks
__global__ void condense_tracks(const device::condense_tracks_payload payload) {
  device::condense_tracks(details::global_index1(), payload);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::condense_tracks_kernel(
    unsigned int n_threads,
    const vecmem::data::vector_view<const unsigned int>&
        out_params_per_in_param,
    vecmem::data::vector_view<unsigned int>& params_index,
    const device::condense_tracks_payload& payload) const {
  // Get a convenience variable for the stream that we'll be using.
  hipStream_t str = details::get_stream(stream());

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

  // Prepare the "in_params_index_view" input parameter of the track
  // condensing kernel, with Thrust's help.
  const vecmem::device_vector<const unsigned int>
      out_params_per_in_param_vector(out_params_per_in_param);
  vecmem::device_vector<unsigned int> params_index_vector(params_index);
  thrust::inclusive_scan(policy, out_params_per_in_param_vector.begin(),
                         out_params_per_in_param_vector.end(),
                         params_index_vector.begin());

  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = warp_size() * 8;
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;

  // Launch the kernel.
  assert(params_index.ptr() == payload.in_params_index_view.ptr());
  hipLaunchKernelGGL(kernels::condense_tracks, deviceBlocks, deviceThreads, 0,
                     str, payload);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

}  // namespace traccc::hip
