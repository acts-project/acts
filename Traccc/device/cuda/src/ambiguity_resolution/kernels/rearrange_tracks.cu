/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// Local include(s).
#include "../../utils/global_index.hpp"
#include "rearrange_tracks.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::cuda::kernels {

TRACCC_DEVICE inline bool find_valid_index(
    unsigned int& idx, const int lower_bound, const int upper_bound,
    const vecmem::device_vector<const unsigned int>& sorted_ids,
    const vecmem::device_vector<const int>& is_updated) {

    const auto initial_idx = idx;

    for (int i = initial_idx; i <= upper_bound; i++) {
        if (!is_updated[sorted_ids[i]]) {
            idx = i;
            return true;
        }
    }

    for (int i = initial_idx - 1; i >= lower_bound; i--) {
        if (!is_updated[sorted_ids[i]]) {
            idx = i;
            return true;
        }
    }

    return false;
}

__launch_bounds__(1024) __global__
    void rearrange_tracks(device::rearrange_tracks_payload payload) {

    if (*(payload.terminate) == 1 || *(payload.n_updated_tracks) == 0) {
        return;
    }

    auto gid = threadIdx.x / nThreads_per_track +
               blockIdx.x * (blockDim.x / nThreads_per_track);
    const unsigned int n_accepted = *(payload.n_accepted);

    auto N = *(payload.n_updated_tracks);

    int neff_threads = (N + nThreads_per_track - 1) / nThreads_per_track;

    if (neff_threads > nThreads_per_track) {
        neff_threads = nThreads_per_track;
    }

    bool is_valid_thread = true;
    if (threadIdx.x % nThreads_per_track >= neff_threads || gid >= n_accepted) {
        is_valid_thread = false;
    }

    vecmem::device_vector<const unsigned int> sorted_ids(
        payload.sorted_ids_view);
    vecmem::device_vector<const unsigned int> inverted_ids(
        payload.inverted_ids_view);
    vecmem::device_vector<const traccc::scalar> rel_shared(
        payload.rel_shared_view);
    vecmem::device_vector<const traccc::scalar> pvals(payload.pvals_view);
    vecmem::device_vector<const unsigned int> updated_tracks(
        payload.updated_tracks_view);
    vecmem::device_vector<const int> is_updated(payload.is_updated_view);
    vecmem::device_vector<const int> prefix_sums(payload.prefix_sums_view);
    vecmem::device_vector<unsigned int> temp_sorted_ids(
        payload.temp_sorted_ids_view);

    __shared__ int shifted_indices[1024];
    auto& shifted_idx = shifted_indices[threadIdx.x / nThreads_per_track];
    unsigned int tid = std::numeric_limits<unsigned int>::max();

    if (is_valid_thread) {

        tid = sorted_ids[gid];
        auto rel_sh_ref = rel_shared[tid];
        auto pval_ref = pvals[tid];

        shifted_idx = static_cast<int>(gid);

        int stride = (N + neff_threads - 1) / neff_threads;

        int ini_idx = stride * (threadIdx.x % nThreads_per_track);
        int fin_idx = std::min(ini_idx + stride, static_cast<int>(N));

        // If it is an updated track, find new sorted index by using the binary
        // search. The index is also shifted by using the bitonic sort result
        // from sort_updated_tracks and prefix sums
        if (is_updated[tid]) {

            if (gid > 0) {

                unsigned int left = 0;
                unsigned int right = gid;

                bool first_iteration = true;

                if (threadIdx.x % nThreads_per_track == 0) {

                    while (right > left) {

                        const bool find_left = find_valid_index(
                            left, 0, gid, sorted_ids, is_updated);

                        if (!find_left) {
                            break;
                        }

                        const bool find_right = find_valid_index(
                            right, 0, gid, sorted_ids, is_updated);

                        if (!find_right) {
                            break;
                        }

                        if (first_iteration) {
                            const auto right_idx = sorted_ids[right];
                            auto rel_sh = rel_shared[right_idx];
                            auto pval = pvals[right_idx];

                            if (rel_sh < rel_sh_ref ||
                                (rel_sh == rel_sh_ref && pval >= pval_ref)) {
                                left = gid;
                                break;
                            }
                        }

                        first_iteration = false;

                        unsigned int mid = left + (right - left) / 2;

                        const bool find_mid = find_valid_index(
                            mid, left, right - 1, sorted_ids, is_updated);

                        if (find_mid) {

                            const auto mid_idx = sorted_ids[mid];
                            auto rel_sh = rel_shared[mid_idx];
                            auto pval = pvals[mid_idx];

                            if (rel_sh < rel_sh_ref ||
                                (rel_sh == rel_sh_ref && pval >= pval_ref)) {

                                left = mid + 1;
                            } else {
                                right = mid;
                            }
                        }
                    }

                    int delta = delta =
                        gid - left - (prefix_sums[gid] - prefix_sums[left]);

                    if (!is_updated[sorted_ids[left]]) {
                        delta++;
                    }

                    atomicAdd(&shifted_idx, -delta);
                }
            }

            for (int i = ini_idx; i < fin_idx; i++) {

                auto id = updated_tracks[i];

                if (inverted_ids[id] < gid) {
                    atomicAdd(&shifted_idx, -1);
                }
            }

            int offset = 0;
            for (int i = ini_idx; i < fin_idx; i++) {
                if (updated_tracks[i] == tid) {
                    offset = i;
                    break;
                }
            }
            if (offset != 0) {
                atomicAdd(&shifted_idx, offset);
            }
        }
        // If it is not an updated track, it is enough to count the number of
        // updated tracks which need to come earlier.
        else {

            for (int i = ini_idx; i < fin_idx; i++) {

                auto id = updated_tracks[i];
                auto rel_sh = rel_shared[id];
                auto pval = pvals[id];

                if (inverted_ids[id] > gid) {
                    if (rel_sh < rel_sh_ref) {
                        atomicAdd(&shifted_idx, 1);
                    } else if (rel_sh == rel_sh_ref && pval > pval_ref) {
                        atomicAdd(&shifted_idx, 1);
                    }
                }
            }
        }
    }

    __syncthreads();

    // Save the result of new indices into a temporary buffer
    if (is_valid_thread && (threadIdx.x % nThreads_per_track) == 0) {
        temp_sorted_ids.at(shifted_idx) = tid;
    }
}

}  // namespace traccc::cuda::kernels
