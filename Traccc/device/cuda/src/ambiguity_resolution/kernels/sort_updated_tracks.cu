/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/barrier.hpp"
#include "../../utils/global_index.hpp"
#include "sort_updated_tracks.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

__launch_bounds__(512) __global__
    void sort_updated_tracks(device::sort_updated_tracks_payload payload) {

    if (*(payload.terminate) == 1 || *(payload.n_updated_tracks) == 0) {
        return;
    }

    __shared__ unsigned int shared_mem_tracks[512];

    vecmem::device_vector<const traccc::scalar> rel_shared(
        payload.rel_shared_view);
    vecmem::device_vector<const traccc::scalar> pvals(payload.pvals_view);
    vecmem::device_vector<unsigned int> updated_tracks(
        payload.updated_tracks_view);

    const unsigned int tid = threadIdx.x;

    // Load updated track indices into shared memory (for sorting)
    shared_mem_tracks[tid] = std::numeric_limits<unsigned int>::max();

    if (tid < *(payload.n_updated_tracks)) {
        shared_mem_tracks[tid] = updated_tracks[tid];
    }

    __syncthreads();

    // Padding the number of tracks to the power of 2
    const unsigned int N = 1 << (32 - __clz(*(payload.n_updated_tracks) - 1));

    traccc::scalar rel_i;
    traccc::scalar rel_j;
    traccc::scalar pval_i;
    traccc::scalar pval_j;

    // Bitonic sort
    for (int k = 2; k <= N; k <<= 1) {

        bool ascending = ((tid & k) == 0);

        for (int j = k >> 1; j > 0; j >>= 1) {
            int ixj = tid ^ j;

            if (ixj > tid && ixj < N && tid < N) {
                unsigned int trk_i = shared_mem_tracks[tid];
                unsigned int trk_j = shared_mem_tracks[ixj];

                if (trk_i == std::numeric_limits<unsigned int>::max()) {
                    rel_i = std::numeric_limits<traccc::scalar>::max();
                    pval_i = 0.f;
                } else {
                    rel_i = rel_shared[trk_i];
                    pval_i = pvals[trk_i];
                }

                if (trk_j == std::numeric_limits<unsigned int>::max()) {
                    rel_j = std::numeric_limits<traccc::scalar>::max();
                    pval_j = 0.f;
                } else {
                    rel_j = rel_shared[trk_j];
                    pval_j = pvals[trk_j];
                }

                bool should_swap =
                    (rel_i > rel_j || (rel_i == rel_j && pval_i < pval_j)) ==
                    ascending;

                if (should_swap) {
                    shared_mem_tracks[tid] = trk_j;
                    shared_mem_tracks[ixj] = trk_i;
                }
            }
            __syncthreads();
        }
    }

    // Write back the sorted result from shared memory to global memory
    if (tid < *(payload.n_updated_tracks)) {
        updated_tracks[tid] = shared_mem_tracks[tid];
    }
}

}  // namespace traccc::cuda::kernels
