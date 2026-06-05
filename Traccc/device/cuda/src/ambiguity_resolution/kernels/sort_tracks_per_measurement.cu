/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "sort_tracks_per_measurement.cuh"

// VecMem include(s).
#include <vecmem/containers/jagged_device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void sort_tracks_per_measurement(
    device::sort_tracks_per_measurement_payload payload) {

    __shared__ unsigned int sh_trk_ids[1024];

    vecmem::jagged_device_vector<unsigned int> tracks_per_measurement(
        payload.tracks_per_measurement_view);

    auto tracks = tracks_per_measurement.at(blockIdx.x);
    const unsigned int tid = threadIdx.x;
    const unsigned int n_tracks = tracks.size();

    sh_trk_ids[tid] = std::numeric_limits<unsigned int>::max();

    if (tid < n_tracks) {
        sh_trk_ids[tid] = tracks[tid];
    }

    // Bitonic sort
    const unsigned int N = 1 << (32 - __clz(n_tracks - 1));
    for (int k = 2; k <= N; k <<= 1) {

        bool ascending = ((tid & k) == 0);

        for (int j = k >> 1; j > 0; j >>= 1) {
            int ixj = tid ^ j;

            if (ixj > tid && ixj < N && tid < N) {
                auto trk_i = sh_trk_ids[tid];
                auto trk_j = sh_trk_ids[ixj];

                bool should_swap = (trk_i > trk_j) == ascending;

                if (should_swap) {
                    sh_trk_ids[tid] = trk_j;
                    sh_trk_ids[ixj] = trk_i;
                }
            }
            __syncthreads();
        }
    }

    if (tid < n_tracks) {
        tracks[tid] = sh_trk_ids[tid];
    }
}
}  // namespace traccc::cuda::kernels
