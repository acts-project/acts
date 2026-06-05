/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "update_status.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void update_status(device::update_status_payload payload) {

    if (*(payload.terminate) == 1) {
        return;
    }

    vecmem::device_vector<const unsigned int> temp_sorted_ids(
        payload.temp_sorted_ids_view);
    vecmem::device_vector<unsigned int> sorted_ids(payload.sorted_ids_view);
    vecmem::device_vector<unsigned int> updated_tracks(
        payload.updated_tracks_view);
    vecmem::device_vector<int> is_updated(payload.is_updated_view);
    vecmem::device_vector<const unsigned int> n_shared(payload.n_shared_view);

    auto globalIndex = threadIdx.x + blockIdx.x * blockDim.x;
    const unsigned int n_accepted = *(payload.n_accepted);
    const unsigned int n_updated = *(payload.n_updated_tracks);

    /***********************
     * Update Max Shared
     ***********************/

    unsigned int shared = 0;

    if (globalIndex < n_accepted) {
        auto tid = sorted_ids[globalIndex];
        shared = n_shared[tid];
    }

    for (int offset = 16; offset > 0; offset >>= 1) {
        unsigned int other_shared =
            __shfl_down_sync(0xffffffff, shared, offset);

        if (other_shared > shared) {
            shared = other_shared;
        }
    }

    if (threadIdx.x == 0) {
        atomicMax(payload.max_shared, shared);
    }

    /***********************
     * Update Sorted Ids
     ***********************/

    if (n_updated == 0) {
        return;
    }

    // Reset is_updated vector
    if (globalIndex < n_updated) {
        is_updated[updated_tracks[globalIndex]] = 0;
    }

    if (globalIndex < n_accepted) {
        auto tid = temp_sorted_ids[globalIndex];
        sorted_ids[globalIndex] = tid;
    }
}

}  // namespace traccc::cuda::kernels
