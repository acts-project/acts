/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/finding/device/propagate_to_next_surface.hpp"
#include "traccc/finding/finding_config.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {

template <typename propagator_t, typename bfield_t>
void propagate_to_next_surface(
    const dim3& grid_size, const dim3& block_size, std::size_t shared_mem_size,
    const cudaStream_t& stream, const finding_config& cfg,
    const typename propagator_t::detector_type::const_view_type& det_data,
    const bfield_t& field_data,
    const device::propagate_to_next_surface_payload& payload);

}  // namespace traccc::cuda
