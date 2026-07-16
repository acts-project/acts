/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../../../utils/global_index.hpp"
#include "../../../utils/hints.hpp"
#include "../propagate_to_next_surface.hpp"

// Project include(s).
#include "traccc/finding/device/propagate_to_next_surface.hpp"

namespace traccc::cuda {
namespace kernels {

template <typename propagator_t, typename bfield_t>
__global__ __launch_bounds__(128) void propagate_to_next_surface(
    const __grid_constant__ finding_config cfg,
    const __grid_constant__
    typename propagator_t::detector_type::const_view_type det_data,
    const __grid_constant__ bfield_t field_data,
    const __grid_constant__ device::propagate_to_next_surface_payload payload) {
    // TODO: Reenable this this additional checks for compilation with the ABI
    // enabled.
    // TRACCC_CUDA_SPILL_TO_SHARED_MEMORY;

    device::propagate_to_next_surface<propagator_t>(
        details::global_index1(), cfg, det_data, field_data, payload);
}

}  // namespace kernels

template <typename propagator_t, typename bfield_t>
void propagate_to_next_surface(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream, const finding_config& cfg,
    const typename propagator_t::detector_type::const_view_type& det_data,
    const bfield_t& field_data,
    const device::propagate_to_next_surface_payload& payload) {

    kernels::propagate_to_next_surface<propagator_t>
        <<<grid_size, block_size, shared_mem_size, stream>>>(
            cfg, det_data, field_data, payload);
}
}  // namespace traccc::cuda
