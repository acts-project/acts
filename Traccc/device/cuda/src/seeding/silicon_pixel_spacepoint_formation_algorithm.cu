/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"

// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/seeding/device/form_spacepoints.hpp"

namespace traccc::cuda {
namespace kernels {

/// Kernel wrapping @c device::form_spacepoints
template <typename detector_t>
__global__ void __launch_bounds__(1024, 1)
    form_spacepoints(typename detector_t::view detector,
                     edm::measurement_collection::const_view measurements,
                     edm::spacepoint_collection::view spacepoints)
    requires(traccc::is_detector_traits<detector_t>)
{
    device::form_spacepoints<detector_t>(details::global_index1(), detector,
                                         measurements, spacepoints);
}

}  // namespace kernels

silicon_pixel_spacepoint_formation_algorithm::
    silicon_pixel_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, vecmem::copy& copy,
        cuda::stream& str, std::unique_ptr<const Logger> logger)
    : device::silicon_pixel_spacepoint_formation_algorithm(mr, copy,
                                                           std::move(logger)),
      cuda::algorithm_base(str) {}

void silicon_pixel_spacepoint_formation_algorithm::form_spacepoints_kernel(
    const form_spacepoints_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_measurements + n_threads - 1) / n_threads;
    detector_buffer_visitor<detector_type_list>(
        payload.detector, [&]<typename detector_traits_t>(
                              const typename detector_traits_t::view& det) {
            kernels::form_spacepoints<detector_traits_t>
                <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(
                    det, payload.measurements, payload.spacepoints);
        });
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
