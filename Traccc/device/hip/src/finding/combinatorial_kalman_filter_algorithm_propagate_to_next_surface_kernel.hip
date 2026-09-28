/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/global_index.hpp"
#include "../utils/hip_error_handling.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/finding/device/propagate_to_next_surface.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::propagate_to_next_surface
template <typename propagator_t, typename bfield_t>
__global__ __launch_bounds__(128) void propagate_to_next_surface(
    const finding_config cfg,
    const typename propagator_t::detector_type* const det_data_ptr,
    const bfield_t field_data,
    const device::propagate_to_next_surface_payload payload) {
  device::propagate_to_next_surface<propagator_t>(
      details::global_index1(), cfg, det_data_ptr, field_data, payload);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::propagate_to_next_surface_kernel(
    unsigned int n_threads, const finding_config& config,
    const detector_buffer& detector, const move_only_any& device_detector,
    const magnetic_field& field,
    const device::propagate_to_next_surface_payload& payload) const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = warp_size() * 4;
  const unsigned int deviceBlocks =
      (n_threads + deviceThreads - 1) / deviceThreads;

  // Launch the kernel for the appropriate detector and magnetic field type.
  detector_buffer_magnetic_field_visitor<detector_type_list,
                                         hip::bfield_type_list<scalar>>(
      detector, field,
      [&]<detray::concepts::detector detector_t, typename bfield_view_t>(
          const detray::detector_view_t<detector_t>&,
          const bfield_view_t& bfield) {
        using detector_device_t = detray::detector_device_t<detector_t>;
        using propagator_t =
            traccc::details::ckf_propagator_t<detector_device_t, bfield_view_t>;

        const vecmem::data::vector_buffer<detector_device_t>&
            device_detector_buffer =
                device_detector
                    .as<vecmem::data::vector_buffer<detector_device_t>>();

        // The kernel name carries a top level comma, which the launch macro
        // can not take, so it is passed through a function pointer.
        const auto kernel =
            &kernels::propagate_to_next_surface<propagator_t, bfield_view_t>;
        hipLaunchKernelGGL(kernel, deviceBlocks, deviceThreads, 0,
                           details::get_stream(stream()), config,
                           device_detector_buffer.ptr(), bfield, payload);
        TRACCC_HIP_ERROR_CHECK(hipGetLastError());
      });
}

}  // namespace traccc::hip
