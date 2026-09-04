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
#include "traccc/hip/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"

// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/seeding/device/count_spacepoints.hpp"
#include "traccc/seeding/device/form_spacepoints.hpp"

// (Roc)Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/scan.h>

// System include(s).
#include <cassert>
#include <memory_resource>

namespace traccc::hip {
namespace kernels {

/// Kernel wrapping @c device::count_spacepoints
__global__ void __launch_bounds__(1024, 1) count_spacepoints(
    edm::measurement_collection::const_view measurements,
    vecmem::data::vector_view<unsigned int> spacepoint_flags) {
  device::count_spacepoints(details::global_index1(), measurements,
                            spacepoint_flags);
}

/// Kernel wrapping @c device::form_spacepoints
template <typename detector_t>
__global__ void __launch_bounds__(1024, 1) form_spacepoints(
    typename detector_t::view detector,
    edm::measurement_collection::const_view measurements,
    vecmem::data::vector_view<const unsigned int> spacepoint_index,
    edm::spacepoint_collection::view spacepoints)
  requires(traccc::is_detector_traits<detector_t>)
{
  device::form_spacepoints<detector_t>(details::global_index1(), detector,
                                       measurements, spacepoint_index,
                                       spacepoints);
}

}  // namespace kernels

silicon_pixel_spacepoint_formation_algorithm::
    silicon_pixel_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, const vecmem::copy& copy,
        const stream_wrapper& str, std::unique_ptr<const Logger> logger)
    : device::silicon_pixel_spacepoint_formation_algorithm(mr, copy,
                                                           std::move(logger)),
      hip::algorithm_base(str) {}

void silicon_pixel_spacepoint_formation_algorithm::count_spacepoints_kernel(
    const count_spacepoints_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks =
      (payload.n_measurements + n_threads - 1) / n_threads;
  hipLaunchKernelGGL(kernels::count_spacepoints, n_blocks, n_threads, 0,
                     details::get_stream(stream()), payload.measurements,
                     payload.spacepoint_flags);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
}

void silicon_pixel_spacepoint_formation_algorithm::scan_spacepoint_flags(
    vecmem::data::vector_view<unsigned int>& spacepoint_flags) const {
  assert(spacepoint_flags.size_ptr() == nullptr);

// Set up the Thrust execution policy. Unfortunately this has to be done in
// a backend specific way. :-(
#if defined(__HIP_PLATFORM_AMD__)
  auto policy =
      thrust::hip::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream()));
#elif defined(__HIP_PLATFORM_NVIDIA__)
  auto policy =
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(mr().main)))
          .on(details::get_stream(stream()));
#else
#error "Unrecognised HIP platform."
#endif

  thrust::inclusive_scan(policy, spacepoint_flags.ptr(),
                         spacepoint_flags.ptr() + spacepoint_flags.capacity(),
                         spacepoint_flags.ptr());
}

void silicon_pixel_spacepoint_formation_algorithm::form_spacepoints_kernel(
    const form_spacepoints_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks =
      (payload.n_measurements + n_threads - 1) / n_threads;
  detector_buffer_visitor<detector_type_list>(
      payload.detector, [&]<typename detector_traits_t>(
                            const typename detector_traits_t::view& det) {
        hipLaunchKernelGGL(
            kernels::form_spacepoints<detector_traits_t>, n_blocks, n_threads,
            0, details::get_stream(stream()), det, payload.measurements,
            payload.spacepoint_index, payload.spacepoints);
      });
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());
  // The base class destroys the prefix sum buffer right after this call, so
  // the kernel must finish before returning.
  stream().synchronize();
}

}  // namespace traccc::hip
