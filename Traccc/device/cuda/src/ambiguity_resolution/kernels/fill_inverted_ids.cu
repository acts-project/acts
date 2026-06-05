/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "fill_inverted_ids.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void fill_inverted_ids(device::fill_inverted_ids_payload payload) {

    if (*(payload.terminate) == 1 || *(payload.n_updated_tracks) == 0) {
        return;
    }

    vecmem::device_vector<const unsigned int> sorted_ids(
        payload.sorted_ids_view);
    vecmem::device_vector<unsigned int> inverted_ids(payload.inverted_ids_view);

    auto globalIndex = threadIdx.x + blockIdx.x * blockDim.x;
    const unsigned int n_accepted = *(payload.n_accepted);

    if (globalIndex >= n_accepted) {
        return;
    }

    // Fill the inverted_ids vector which converts a track id to the index of
    // sorted ids
    inverted_ids[sorted_ids[globalIndex]] = globalIndex;
}

}  // namespace traccc::cuda::kernels
