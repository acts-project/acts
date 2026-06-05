/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/barrier.hpp"
#include "../../utils/global_index.hpp"
#include "gather_best_tips_per_measurement.cuh"

namespace traccc::cuda::kernels {

__global__ void gather_best_tips_per_measurement(
    const __grid_constant__
        device::gather_best_tips_per_measurement_payload<default_algebra>
            payload) {

    device::gather_best_tips_per_measurement(details::global_index1(),
                                             barrier{}, payload);
}

}  // namespace traccc::cuda::kernels
