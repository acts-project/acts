/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../../../utils/global_index.hpp"
#include "../propagate_to_next_surface.hpp"

// Project include(s).
#include "traccc/finding/device/propagate_to_next_surface.hpp"

namespace traccc::cuda {
namespace kernels {

template <typename propagator_t, typename bfield_t>
__global__ __launch_bounds__(128) void propagate_to_next_surface(
    const __grid_constant__ finding_config cfg,
    const __grid_constant__
        device::propagate_to_next_surface_payload<propagator_t, bfield_t>
            payload) {

    device::propagate_to_next_surface<propagator_t, bfield_t>(
        details::global_index1(), cfg, payload);
}

}  // namespace kernels

template <typename propagator_t, typename bfield_t>
void propagate_to_next_surface(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream, const finding_config cfg,
    device::propagate_to_next_surface_payload<propagator_t, bfield_t> payload) {

    kernels::propagate_to_next_surface<<<grid_size, block_size, shared_mem_size,
                                         stream>>>(cfg, payload);
}
}  // namespace traccc::cuda
