/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../../../utils/global_index.hpp"
#include "../../../utils/hints.hpp"
#include "../progressive_kalman_filter.hpp"

// Project include(s).
#include "traccc/finding/device/progressive_kalman_filter.hpp"

namespace traccc::cuda {
namespace kernels {

template <typename propagator_t>
__global__ __launch_bounds__(128) void progressive_kalman_filter(
    const __grid_constant__ finding_config cfg,
    const __grid_constant__
    typename propagator_t::detector_type::const_view_type det_data,
    const __grid_constant__
    typename propagator_t::stepper_type::magnetic_field_type field_data,
    const __grid_constant__ vecmem::data::jagged_vector_view<
        typename propagator_t::detector_type::surface_type>
        surfaces_view,
    const __grid_constant__ device::progressive_kalman_filter_payload payload) {
  device::progressive_kalman_filter<propagator_t>(details::global_index1(), cfg,
                                                  det_data, field_data,
                                                  surfaces_view, payload);
}

}  // namespace kernels

template <typename propagator_t>
void progressive_kalman_filter(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream, const finding_config& cfg,
    const typename propagator_t::detector_type::const_view_type& det_data,
    const typename propagator_t::stepper_type::magnetic_field_type& field_data,
    const vecmem::data::jagged_vector_view<
        typename propagator_t::detector_type::surface_type>& surfaces_view,
    const device::progressive_kalman_filter_payload& payload) {
  kernels::progressive_kalman_filter<propagator_t>
      <<<grid_size, block_size, shared_mem_size, stream>>>(
          cfg, det_data, field_data, surfaces_view, payload);
}
}  // namespace traccc::cuda
