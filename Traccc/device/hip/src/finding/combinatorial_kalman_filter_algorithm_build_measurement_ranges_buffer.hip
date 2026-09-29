/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/finding/device/geo_id_surface_comparator.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"

// (Roc)Thrust include(s).
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <memory_resource>

namespace traccc::hip {

vecmem::data::vector_buffer<edm::measurement_collection::const_view::size_type>
combinatorial_kalman_filter_algorithm::build_measurement_ranges_buffer(
    const detector_buffer& det,
    const edm::measurement_collection::const_view::size_type n_measurements,
    const edm::measurement_collection::const_view& measurements) const {
  return detector_buffer_visitor<detector_type_list>(
      det, [&]<detray::concepts::detector detector_t>(
               const detray::detector_view_t<detector_t>& det) {
        // Construct an appropriate device detector object.
        detray::detector_device_t<detector_t> device_det{det};

        // Create the result buffer.
        vecmem::data::vector_buffer<
            edm::measurement_collection::const_view::size_type>
            result{device_det.surfaces().size(), mr().main};
        copy().setup(result)->ignore();

        // Create a measurement device object for convenience.
        const edm::measurement_collection::const_device measurements_device{
            measurements};

        // Get a convenience variable for the stream that we'll be using.
        hipStream_t str = details::get_stream(stream());

// Set up the Thrust execution policy. Unfortunately this has to be done in
// a backend specify way. :-(
#if defined(__HIP_PLATFORM_AMD__)
        auto policy = thrust::hip::par_nosync(
                          std::pmr::polymorphic_allocator(&(mr().main)))
                          .on(str);
#elif defined(__HIP_PLATFORM_NVIDIA__)
        auto policy =
            thrust::cuda::par_nosync(
                std::pmr::polymorphic_allocator(&(mr().main)))
                .on(str);
#else
#error "Unrecognised HIP platform."
#endif

        // Fill it with Thrust's help.
        thrust::upper_bound(
            policy, measurements_device.surface_link().begin(),
            // We have to use this ugly form here, because if the
            // measurement collection is resizable (which it often
            // is), the end() function cannot be used in host code.
            measurements_device.surface_link().begin() + n_measurements,
            device_det.surfaces().begin(), device_det.surfaces().end(),
            result.ptr(), device::geo_id_surface_comparator{});

        // Return the filled buffer.
        return result;
      });
}

}  // namespace traccc::hip
