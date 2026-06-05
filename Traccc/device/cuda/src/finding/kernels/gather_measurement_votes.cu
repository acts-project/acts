/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "gather_measurement_votes.cuh"

namespace traccc::cuda::kernels {

__global__ void gather_measurement_votes(
    const __grid_constant__ device::gather_measurement_votes_payload payload) {

    device::gather_measurement_votes(details::global_index1(), payload);
}

}  // namespace traccc::cuda::kernels
