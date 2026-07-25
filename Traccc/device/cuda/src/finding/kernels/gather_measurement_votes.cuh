/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/finding/device/gather_measurement_votes.hpp"

namespace traccc::cuda::kernels {

__global__ void gather_measurement_votes(
    const __grid_constant__ device::gather_measurement_votes_payload payload);

}  // namespace traccc::cuda::kernels
