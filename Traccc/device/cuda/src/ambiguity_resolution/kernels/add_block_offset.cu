/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/barrier.hpp"
#include "../../utils/global_index.hpp"
#include "add_block_offset.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void add_block_offset(device::add_block_offset_payload payload) {

    if (*(payload.terminate) == 1 || *(payload.n_updated_tracks) == 0) {
        return;
    }

    vecmem::device_vector<const int> block_offsets(payload.block_offsets_view);
    vecmem::device_vector<int> prefix_sums(payload.prefix_sums_view);

    auto globalIndex = threadIdx.x + blockIdx.x * blockDim.x;
    const unsigned int n_accepted = *(payload.n_accepted);

    if (globalIndex >= n_accepted || blockIdx.x == 0) {
        return;
    }

    // Add the scanned block offsets to block-wise prefix sums of the number of
    // updated tracks.
    prefix_sums[globalIndex] += block_offsets[blockIdx.x - 1];
}

}  // namespace traccc::cuda::kernels
