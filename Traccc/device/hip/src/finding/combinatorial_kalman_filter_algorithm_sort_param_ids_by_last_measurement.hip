/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// (Roc)Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/sort.h>

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <cassert>
#include <memory_resource>

namespace traccc::hip {

void combinatorial_kalman_filter_algorithm::sort_param_ids_by_last_measurement(
    vecmem::data::vector_view<unsigned int>& link_last_measurement,
    vecmem::data::vector_view<unsigned int>& param_ids) const {
  assert(link_last_measurement.capacity() == param_ids.capacity());
  assert(link_last_measurement.size_ptr() == nullptr);
  assert(param_ids.size_ptr() == nullptr);

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

  thrust::sort_by_key(
      policy, link_last_measurement.ptr(),
      link_last_measurement.ptr() + link_last_measurement.capacity(),
      param_ids.ptr());
}

}  // namespace traccc::hip
