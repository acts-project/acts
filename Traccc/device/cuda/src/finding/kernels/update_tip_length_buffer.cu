/** traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "update_tip_length_buffer.cuh"

namespace traccc::cuda::kernels {

__global__ void update_tip_length_buffer(
    const __grid_constant__ device::update_tip_length_buffer_payload payload) {

    device::update_tip_length_buffer(details::global_index1(), payload);
}

}  // namespace traccc::cuda::kernels
