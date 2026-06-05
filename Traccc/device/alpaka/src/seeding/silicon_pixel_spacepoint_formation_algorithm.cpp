/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/alpaka/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"

#include "../utils/barrier.hpp"
#include "../utils/get_queue.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/seeding/device/form_spacepoints.hpp"

namespace traccc::alpaka {
namespace kernels {

/// Kernel for running @c traccc::device::form_spacepoints
template <typename detector_t>
struct form_spacepoints {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const typename detector_t::view* detector,
        edm::measurement_collection::const_view measurements,
        edm::spacepoint_collection::view spacepoints) const {

        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];

        device::form_spacepoints<detector_t>(globalThreadIdx, *detector,
                                             measurements, spacepoints);
    }
};

}  // namespace kernels

silicon_pixel_spacepoint_formation_algorithm::
    silicon_pixel_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, vecmem::copy& copy, alpaka::queue& q,
        std::unique_ptr<const Logger> logger)
    : device::silicon_pixel_spacepoint_formation_algorithm(mr, copy,
                                                           std::move(logger)),
      alpaka::algorithm_base(q) {}

void silicon_pixel_spacepoint_formation_algorithm::form_spacepoints_kernel(
    const form_spacepoints_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 8;
    const unsigned int n_blocks =
        (payload.n_measurements + n_threads - 1) / n_threads;
    detector_buffer_visitor<detector_type_list>(
        payload.detector, [&]<typename detector_traits_t>(
                              const typename detector_traits_t::view& det) {
            // Put the detector view into device memory.
            vecmem::data::vector_buffer<typename detector_traits_t::view>
                device_det(1u, mr().main);
            copy().setup(device_det)->wait();
            copy()(vecmem::data::vector_view<
                       const typename detector_traits_t::view>(1u, &det),
                   device_det)
                ->wait();
            // Launch the spacepoint formation kernel.
            ::alpaka::exec<Acc>(details::get_queue(queue()),
                                makeWorkDiv<Acc>(n_blocks, n_threads),
                                kernels::form_spacepoints<detector_traits_t>{},
                                device_det.ptr(), payload.measurements,
                                payload.spacepoints);
        });
}

}  // namespace traccc::alpaka
