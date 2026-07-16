/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/barrier.hpp"
#include "../../utils/global_index.hpp"
#include "scan_block_offsets.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void scan_block_offsets(device::scan_block_offsets_payload payload) {

    if (*(payload.terminate) == 1 || *(payload.n_updated_tracks) == 0) {
        return;
    }

    extern __shared__ int shared_temp[];

    vecmem::device_vector<const int> block_offsets(payload.block_offsets_view);
    vecmem::device_vector<int> scanned_block_offsets(
        payload.scanned_block_offsets_view);

    // The number of blocks in the previous block_inclusive_scan = the nubmer of
    // threads of this kernel
    int n_blocks_prev = blockDim.x;
    auto threadIndex = threadIdx.x;

    // 1. Load from global to shared
    shared_temp[threadIndex] = 0;
    if (threadIndex < n_blocks_prev) {
        shared_temp[threadIndex] = block_offsets[threadIndex];
    }
    __syncthreads();

    // 2. Inclusive scan to caculated the scanned block offset which is the
    // prefix sum of block offsets
    for (int offset = 1; offset < n_blocks_prev; offset *= 2) {
        int temp = 0;
        if (threadIndex >= offset) {
            temp = shared_temp[threadIndex - offset];
        }
        __syncthreads();
        shared_temp[threadIndex] += temp;
        __syncthreads();
    }

    // 3. Write back
    if (threadIndex < n_blocks_prev) {
        scanned_block_offsets[threadIndex] = shared_temp[threadIndex];
    }
}

}  // namespace traccc::cuda::kernels
