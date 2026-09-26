/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
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
#include "traccc/finding/details/progressive_kalman_filter_types.hpp"
#include "traccc/finding/device/progressive_kalman_filter.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"
#include "traccc/utils/propagation.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {
namespace kernels {

/// HIP kernel for running @c traccc::device::progressive_kalman_filter
template <typename propagator_t>
__global__ __launch_bounds__(128) void progressive_kalman_filter(
    const finding_config cfg,
    const typename propagator_t::detector_type::const_view_type det_data,
    const typename propagator_t::stepper_type::magnetic_field_type field_data,
    const vecmem::data::jagged_vector_view<
        typename propagator_t::detector_type::surface_type>
        surfaces_view,
    const device::progressive_kalman_filter_payload payload) {
  device::progressive_kalman_filter<propagator_t>(details::global_index1(), cfg,
                                                  det_data, field_data,
                                                  surfaces_view, payload);
}

}  // namespace kernels

void combinatorial_kalman_filter_algorithm::progressive_kalman_filter_kernel(
    unsigned int n_seeds, const finding_config& config,
    const detector_buffer& detector, const magnetic_field& field,
    const device::progressive_kalman_filter_payload& payload,
    const device::kalman_fitting_algorithm::fit_payload& smoothing_payload)
    const {
  // Establish the kernel launch parameters.
  const unsigned int deviceThreads = 128;
  const unsigned int deviceBlocks =
      (n_seeds + deviceThreads - 1) / deviceThreads;

  detector_buffer_magnetic_field_visitor<detector_type_list,
                                         hip::bfield_type_list<scalar>>(
      detector, field,
      [&]<detray::concepts::detector detector_t, typename bfield_view_t>(
          const detray::detector_view_t<detector_t>& det,
          const bfield_view_t& bfield) {
        using detector_device_t = detray::detector_device_t<detector_t>;
        using surface_t = typename detector_device_t::surface_type;
        using propagator_t =
            traccc::details::pkf_propagator_t<detector_device_t, bfield_view_t>;

        // If the Kalman smoother should be run, obtain the real allocation
        vecmem::data::jagged_vector_view<surface_t> sf_sequences;
        if (config.run_smoother == smoother_type::e_kalman) {
          sf_sequences =
              smoothing_payload.surfaces
                  .as<vecmem::data::jagged_vector_buffer<surface_t>>();
        }

        // Launch the kernel.
        hipLaunchKernelGGL(kernels::progressive_kalman_filter<propagator_t>,
                           deviceBlocks, deviceThreads, 0,
                           details::get_stream(stream()), config, det, bfield,
                           sf_sequences, payload);
        TRACCC_HIP_ERROR_CHECK(hipGetLastError());
      });
}

}  // namespace traccc::hip
