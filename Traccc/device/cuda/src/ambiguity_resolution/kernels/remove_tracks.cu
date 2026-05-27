/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/utils/pair.hpp"

// Local include(s).
#include "../../utils/barrier.hpp"
#include "../../utils/global_index.hpp"
#include "remove_tracks.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

// Thrust include(s).
#include <thrust/binary_search.h>
#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/find.h>

namespace traccc::cuda::kernels {

__device__ __forceinline__ uint64_t pack_key(uint32_t meas, uint32_t thr) {
    return (uint64_t(meas) << 32) | uint64_t(thr);
}
__device__ __forceinline__ void unpack_key(uint64_t k, uint32_t& meas,
                                           uint32_t& thr) {
    meas = uint32_t(k >> 32);
    thr = uint32_t(k & 0xFFFFFFFFu);
}

__device__ void count_tracks(int tid, int* sh_n_meas, int n_tracks,
                             unsigned int& bound, unsigned int& count) {

    unsigned int add = 0;

    // --- Warp-level phase: handle strides < 32 using warp shuffle (no
    // __syncthreads needed) ---
    const int lane = threadIdx.x & 31;
    const unsigned int full_mask = 0xFFFFFFFFu;

    // Load this thread's value into a register if it's in range
    int v = (tid < n_tracks) ? sh_n_meas[tid] : 0;

    // Mask for active lanes in this warp
    const unsigned int mask = __ballot_sync(full_mask, tid < n_tracks);

    const int max_stride = min(n_tracks, 32);

    for (int stride = 1; stride < max_stride; stride <<= 1) {
        // Accumulate neighbor's value via warp shuffle
        unsigned int other = __shfl_down_sync(mask, v, stride);
        if (lane + stride < 32 && (tid + stride) < n_tracks) {
            v += other;
        }

        // Thread 0 can directly check its register value in the warp phase
        if (tid == 0) {
            if (v < bound) {
                add = stride << 1;
            }
        }
    }

    // Write warp-phase result back to shared memory
    if (tid < n_tracks) {
        sh_n_meas[tid] = static_cast<int>(v);
    }
    __syncthreads();

    // --- Block-level phase: handle strides >= 32 (minimal required
    // synchronizations) ---
    for (int stride = 32; stride < n_tracks; stride <<= 1) {
        if ((tid + stride) < n_tracks) {
            sh_n_meas[tid] += sh_n_meas[tid + stride];
        }
        __syncthreads();

        if (tid == 0 && sh_n_meas[0] < bound) {
            add = stride << 1;
        }
        __syncthreads();
    }

    // --- Final update ---
    if (tid == 0) {
        count += add;
    }
    __syncthreads();
}

__launch_bounds__(512) __global__
    void remove_tracks(device::remove_tracks_payload payload) {

    if (threadIdx.x == 0) {
        if (*(payload.max_shared) == 0) {
            *(payload.terminate) = 1;
        }
        // Reset the max_shared and n_updated_tracks
        *(payload.max_shared) = 0;
        *(payload.n_updated_tracks) = 0;
    }

    __syncthreads();

    if (*(payload.terminate) == 1) {
        return;
    }

    __shared__ int sh_buffer[512];
    __shared__ measurement_id_type sh_meas_ids[512];
    __shared__ unsigned int sh_threads[512];
    __shared__ uint64_t sh_keys[512];
    __shared__ unsigned int n_meas_total;
    __shared__ unsigned int bound;
    __shared__ unsigned int n_tracks_to_iterate;
    __shared__ unsigned int min_thread;
    __shared__ unsigned int N;
    __shared__ unsigned int n_updating_threads;

    auto threadIndex = threadIdx.x;

    int gid = static_cast<int>(*payload.n_accepted) - 1 - threadIndex;
    sh_buffer[threadIndex] = 0;
    sh_meas_ids[threadIndex] = std::numeric_limits<measurement_id_type>::max();
    sh_threads[threadIndex] = std::numeric_limits<unsigned int>::max();

    vecmem::device_vector<const unsigned int> sorted_ids(
        payload.sorted_ids_view);
    vecmem::jagged_device_vector<const measurement_id_type> meas_ids(
        payload.meas_ids_view);
    vecmem::device_vector<const unsigned int> n_meas(payload.n_meas_view);
    vecmem::device_vector<const unsigned int> meas_id_to_unique_id(
        payload.meas_id_to_unique_id_view);
    vecmem::jagged_device_vector<const unsigned int> tracks_per_measurement(
        payload.tracks_per_measurement_view);
    vecmem::jagged_device_vector<int> track_status_per_measurement(
        payload.track_status_per_measurement_view);
    vecmem::device_vector<unsigned int> n_accepted_tracks_per_measurement(
        payload.n_accepted_tracks_per_measurement_view);
    vecmem::device_vector<unsigned int> n_shared(payload.n_shared_view);
    vecmem::device_vector<traccc::scalar> rel_shared(payload.rel_shared_view);
    vecmem::device_vector<unsigned int> updated_tracks(
        payload.updated_tracks_view);
    vecmem::device_vector<int> is_updated(payload.is_updated_view);
    vecmem::device_vector<int> track_count(payload.track_count_view);

    if (threadIndex == 0) {
        *(payload.n_removable_tracks) = 0;
        *(payload.n_meas_to_remove) = 0;
        *(payload.n_valid_threads) = 0;
        n_meas_total = 0;
        bound = 512;
        N = 1;
        n_tracks_to_iterate = 0;
        min_thread = std::numeric_limits<unsigned int>::max();
    }

    __syncthreads();

    unsigned int trk_id = 0;
    unsigned int n_m = 0;
    if (gid >= 0) {
        trk_id = sorted_ids[gid];
        n_m = n_meas[trk_id];

        // Buffer for the number of measurement per track
        sh_buffer[threadIndex] = n_m;
    }

    __syncthreads();

    auto n_tracks_total = min(bound, *payload.n_accepted);

    /****************************************
     * Count the number of removable tracks
     ****************************************/

    count_tracks(threadIdx.x, sh_buffer, n_tracks_total, bound,
                 n_tracks_to_iterate);

    if (threadIndex == 0 && n_tracks_to_iterate == 0) {
        n_tracks_to_iterate = 1;
    }

    // @TODO: Improve the logic
    if (threadIndex < n_tracks_to_iterate && gid >= 0) {
        const unsigned int pos = atomicAdd(&n_meas_total, n_m);

        const auto& mids = meas_ids[trk_id];
        for (int i = 0; i < n_m; i++) {
            sh_meas_ids[pos + i] = mids[i];
            sh_threads[pos + i] = threadIndex;
        }
    }

    __syncthreads();

    // Bitonic sort on meas_to_thread w.r.t. measurement id
    if (threadIndex == 0) {
        N = (n_meas_total == 0) ? 1 : 1 << (32 - __clz(n_meas_total - 1));
    }
    __syncthreads();

    const auto tid = threadIndex;
    // No early return: out-of-range threads carry a sentinel and only
    // sync/shuffle.
    uint64_t key = (tid < N) ? pack_key(sh_meas_ids[tid], sh_threads[tid])
                             : 0xFFFFFFFFFFFFFFFFull;  // sentinel that won't
                                                       // affect in-range items

    for (int k = 2; k <= N; k <<= 1) {
        // Inter-warp (j >= 32): use shared + barriers
        for (int j = (k >> 1); j >= warpSize; j >>= 1) {
            sh_keys[tid] = key;  // safe: sh_keys sized to blockDim.x
            __syncthreads();

            const int ixj = tid ^ j;
            // If partner is out-of-range, compare with self (no change).
            uint64_t other = (ixj < N) ? sh_keys[ixj] : key;

            const bool dir = ((tid & k) == 0);    // ascending segment?
            const bool lower = ((tid & j) == 0);  // am I lower index?

            const uint64_t mn = (key < other) ? key : other;
            const uint64_t mx = (key < other) ? other : key;

            key = dir ? (lower ? mn : mx) : (lower ? mx : mn);

            __syncthreads();
        }

        // Intra-warp (j < 32): warp shuffles only; no barriers
        for (int j = min(k >> 1, warpSize >> 1); j > 0; j >>= 1) {
            const unsigned mask = 0xFFFFFFFFu;
            uint64_t other = __shfl_xor_sync(mask, key, j);

            const bool dir = ((tid & k) == 0);
            const bool lower = ((tid & j) == 0);

            const uint64_t mn = (key < other) ? key : other;
            const uint64_t mx = (key < other) ? other : key;

            key = dir ? (lower ? mn : mx) : (lower ? mx : mn);
        }

        // Commit for next inter-warp round visibility
        sh_keys[tid] = key;
        __syncthreads();
    }

    // Write back only in-range threads
    if (tid < N) {
        uint32_t meas, thr;
        unpack_key(key, meas, thr);
        sh_meas_ids[tid] = meas;
        sh_threads[tid] = thr;
    }

    // Find starting point
    if (threadIndex < n_meas_total) {
        auto mid = sh_meas_ids[threadIndex];
        bool is_start =
            (threadIndex == 0) || (sh_meas_ids[threadIndex - 1] != mid);
        const auto unique_meas_idx = meas_id_to_unique_id.at(mid);
        const auto its_accepted_tracks =
            n_accepted_tracks_per_measurement.at(unique_meas_idx);

        if (is_start) {

            int i = threadIndex + 1;
            int n_sharing_tracks = 1;

            while (i < n_meas_total && sh_meas_ids[i] == mid) {
                if (sh_threads[i] != sh_threads[i - 1]) {
                    n_sharing_tracks++;

                    if (n_sharing_tracks == its_accepted_tracks) {
                        atomicMin(&min_thread, sh_threads[i - 1]);
                        break;
                    }
                }
                i++;
            }
        }
    }

    __syncthreads();

    if (threadIndex == 0) {
        if (min_thread == 0) {
            *(payload.n_removable_tracks) = 1;
        } else if (min_thread == std::numeric_limits<unsigned int>::max()) {
            *(payload.n_removable_tracks) = n_tracks_to_iterate;
        } else {
            *(payload.n_removable_tracks) = min_thread;
        }
    }

    __syncthreads();

    auto meas_to_remove_temp = sh_meas_ids[threadIndex];
    auto threads_temp = sh_threads[threadIndex];

    __syncthreads();

    if (threadIndex == 0) {
        *(payload.n_meas_to_remove) = n_meas_total;
    }

    __syncthreads();

    int is_valid = (threads_temp < *(payload.n_removable_tracks)) ? 1 : 0;

    // TODO: Use better reduction algorithm
    if (is_valid) {
        atomicAdd(payload.n_valid_threads, 1);
    }

    __syncthreads();

    // Exclusive scan (Hillis-Steele)

    // Buffer for the prefix
    sh_buffer[threadIndex] = is_valid;  // copy input
    __syncthreads();

    for (int offset = 1; offset < *(payload.n_meas_to_remove); offset <<= 1) {
        int val = 0;
        if (threadIndex >= offset) {
            val = sh_buffer[threadIndex - offset];
        }
        __syncthreads();
        sh_buffer[threadIndex] += val;
        __syncthreads();
    }

    if (is_valid) {
        sh_buffer[threadIndex] -= 1;
        sh_meas_ids[sh_buffer[threadIndex]] = meas_to_remove_temp;
        sh_threads[sh_buffer[threadIndex]] = threads_temp;
    }

    __syncthreads();

    meas_to_remove_temp = sh_meas_ids[threadIndex];
    threads_temp = sh_threads[threadIndex];

    /********************
     * Remove tracks
     ********************/

    __syncthreads();

    bool is_valid_thread = false;
    bool is_duplicate = true;

    const unsigned n_accepted_prev = *(payload.n_accepted);

    __syncthreads();

    if (threadIndex == 0) {
        (*payload.n_accepted) -= *(payload.n_removable_tracks);
        n_updating_threads = 0;
    }

    if (threadIndex < *(payload.n_valid_threads)) {
        sh_meas_ids[threadIndex] = meas_to_remove_temp;
        sh_threads[threadIndex] = threads_temp;
        is_valid_thread = true;
    }

    __syncthreads();

    if (is_valid_thread) {
        const auto id = sh_meas_ids[threadIndex];
        is_duplicate = (threadIndex > 0 && sh_meas_ids[threadIndex - 1] == id);
    }

    // Buffer for the track ids
    sh_buffer[threadIndex] = std::numeric_limits<int>::max();

    bool active = false;
    unsigned int pos1;
    int alive_trk_id = 0;

    if (!is_duplicate && is_valid_thread) {

        const auto id = sh_meas_ids[threadIndex];
        const auto unique_meas_idx = meas_id_to_unique_id.at(id);

        // If there is only one track associated with measurement, the
        // number of shared measurement can be reduced by one
        const auto& tracks = tracks_per_measurement[unique_meas_idx];
        auto track_status = track_status_per_measurement[unique_meas_idx];

        auto trk_id =
            sorted_ids.at(n_accepted_prev - 1 - sh_threads[threadIndex]);

        unsigned int worst_idx =
            thrust::lower_bound(thrust::seq, tracks.begin(), tracks.end(),
                                trk_id) -
            tracks.begin();

        track_status[worst_idx] = 0;

        int n_sharing_tracks = 1;
        for (unsigned int i = threadIndex + 1; i < *(payload.n_valid_threads);
             ++i) {

            if (sh_meas_ids[i] == id && sh_threads[i] != sh_threads[i - 1]) {
                n_sharing_tracks++;

                trk_id = sorted_ids[n_accepted_prev - 1 - sh_threads[i]];

                worst_idx = thrust::lower_bound(thrust::seq, tracks.begin(),
                                                tracks.end(), trk_id) -
                            tracks.begin();

                track_status[worst_idx] = 0;

            } else if (sh_meas_ids[i] != id) {
                break;
            }
        }

        vecmem::device_atomic_ref<unsigned int> n_accepted_per_meas(
            n_accepted_tracks_per_measurement.at(
                static_cast<unsigned int>(unique_meas_idx)));
        const unsigned int N_A =
            n_accepted_per_meas.fetch_sub(n_sharing_tracks);

        if (N_A == 1 + n_sharing_tracks) {
            active = true;
            const unsigned int alive_idx =
                thrust::find(thrust::seq, track_status.begin(),
                             track_status.end(), 1) -
                track_status.begin();

            pos1 = atomicAdd(&n_updating_threads, 1);
            alive_trk_id = static_cast<int>(tracks[alive_idx]);

            sh_buffer[pos1] = alive_trk_id;
            atomicAdd(&track_count[alive_trk_id], 1);

            const auto m_count = static_cast<unsigned int>(
                thrust::count(thrust::seq, meas_ids[alive_trk_id].begin(),
                              meas_ids[alive_trk_id].end(), id));

            const unsigned int N_S = vecmem::device_atomic_ref<unsigned int>(
                                         n_shared.at(alive_trk_id))
                                         .fetch_sub(m_count);
        }
    }

    __syncthreads();

    if (active) {

        auto count = atomicAdd(&track_count[alive_trk_id], -1);
        if (count == 1) {

            // Write updated track IDs
            vecmem::device_atomic_ref<unsigned int> num_updated_tracks(
                *(payload.n_updated_tracks));

            const unsigned int pos2 = num_updated_tracks.fetch_add(1);

            updated_tracks[pos2] = alive_trk_id;
            is_updated[alive_trk_id] = 1;

            rel_shared.at(alive_trk_id) = math::div_ieee754(
                static_cast<traccc::scalar>(n_shared.at(alive_trk_id)),
                static_cast<traccc::scalar>(n_meas.at(alive_trk_id)));
        }
    }
}

}  // namespace traccc::cuda::kernels
