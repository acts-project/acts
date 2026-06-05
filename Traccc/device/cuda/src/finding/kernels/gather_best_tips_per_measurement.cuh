/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/finding/device/gather_best_tips_per_measurement.hpp"

namespace traccc::cuda::kernels {

__global__ void gather_best_tips_per_measurement(
    const __grid_constant__
        device::gather_best_tips_per_measurement_payload<default_algebra>
            payload);

}  // namespace traccc::cuda::kernels
