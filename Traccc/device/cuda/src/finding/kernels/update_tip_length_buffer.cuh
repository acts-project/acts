/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/finding/device/update_tip_length_buffer.hpp"

namespace traccc::cuda::kernels {

__global__ void update_tip_length_buffer(
    const __grid_constant__ device::update_tip_length_buffer_payload payload);

}  // namespace traccc::cuda::kernels
