/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "../../utils/global_index.hpp"
#include "remove_duplicates.cuh"
#include "traccc/finding/device/remove_duplicates.hpp"

namespace traccc::cuda::kernels {

__global__ void remove_duplicates(
    const __grid_constant__ finding_config cfg,
    const __grid_constant__ device::remove_duplicates_payload payload) {

    device::remove_duplicates(details::global_index1(), cfg, payload);
}

}  // namespace traccc::cuda::kernels
