/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/barrier.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t, concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void gbts_fill_path_store(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const gbts_fill_path_store_payload& payload,
    const gbts_fill_path_store_shared_payload& shared) {

    const unsigned int threadIndex = thread_id.getLocalThreadIdX();
    const unsigned int blockIndex = thread_id.getBlockIdX();
    const unsigned int blockSize = thread_id.getBlockDimX();

    vecmem::device_vector<int2> d_path_store(payload.path_store);
    const vecmem::device_vector<const unsigned int> d_output_graph(
        payload.output_graph);
    const vecmem::device_vector<const unsigned char> d_levels(payload.levels);

    vecmem::device_vector<traccc::uint2> shared_live_paths(shared.live_paths);

    unsigned int& nPathStoreSizeCounter = *payload.nPathStoreSizeCounter;

    // Number of terminus edges expanded per block: bounded by the block size
    // and by how many of the (on average nPaths / nTerminusEdges) paths per
    // terminus fit into the shared live-path buffer.
    const unsigned int pathsPerTerminus =
        1u + (payload.nPaths - 1u) / payload.nTerminusEdges;
    const unsigned int bufferTerminusPerBlock =
        1u + (static_cast<unsigned int>(
                  traccc::device::gbts_consts::live_path_buffer) -
              1u) /
                 pathsPerTerminus;
    const unsigned int nTerminusPerBlock = (blockSize < bufferTerminusPerBlock)
                                               ? blockSize
                                               : bufferTerminusPerBlock;
    // Upper bound on path-store slots that may be claimed.
    const unsigned int nReachablePaths =
        payload.nPaths + payload.nTerminusEdges;

    if (threadIndex == 0) {
        shared.n_live_paths = 0;
    }
    barrier.blockBarrier();

    // Row-major output graph: each edge owns a contiguous block of
    // edge_size = 2 + 1 + max_num_neighbours ints.
    const unsigned int edge_size = 2u + 1u + payload.max_num_neighbours;
    unsigned int path_idx = threadIndex + blockIndex * nTerminusPerBlock;

    if (threadIndex < nTerminusPerBlock && path_idx < payload.nTerminusEdges) {
        const int2 path = d_path_store[path_idx];
        const unsigned int edge_pos =
            edge_size * static_cast<unsigned int>(path.x);
        const unsigned int nNei = d_output_graph[edge_pos + gbts_consts::nNei];
        const unsigned char level = d_levels[static_cast<unsigned int>(path.x)];
        for (unsigned int nei = 0; nei < nNei; ++nei) {
            const unsigned int edge_idx =
                d_output_graph[edge_pos + gbts_consts::nei_start + nei];
            if (level != d_levels[static_cast<unsigned int>(edge_idx)] + 1) {
                continue;
            }
            const unsigned int live_idx = static_cast<unsigned int>(
                vecmem::device_atomic_ref<int,
                                          vecmem::device_address_space::local>(
                    shared.n_live_paths)
                    .fetch_add(1));
            if (live_idx >=
                static_cast<unsigned int>(
                    traccc::device::gbts_consts::live_path_buffer)) {
                break;
            }
            const unsigned int new_path_idx =
                vecmem::device_atomic_ref<unsigned int>(nPathStoreSizeCounter)
                    .fetch_add(1u);
            d_path_store[new_path_idx] =
                int2{static_cast<int>(edge_idx), static_cast<int>(path_idx)};
            shared_live_paths[static_cast<unsigned int>(live_idx)] =
                uint2{edge_idx, new_path_idx};
        }
    }
    barrier.blockBarrier();

    traccc::uint2 path = uint2{0u, 0u};
    bool has_path = false;

    while (shared.n_live_paths > 0) {
        has_path = false;
        if (threadIndex == 0) {
            const int buf_size =
                static_cast<int>(traccc::device::gbts_consts::live_path_buffer);
            shared.n_live_paths = (shared.n_live_paths < buf_size)
                                      ? shared.n_live_paths
                                      : buf_size;
        }
        barrier.blockBarrier();
        if (static_cast<int>(threadIndex) < shared.n_live_paths) {
            path = shared_live_paths[static_cast<unsigned int>(
                shared.n_live_paths - static_cast<int>(threadIndex) - 1)];
            has_path = true;
        }
        barrier.blockBarrier();
        if (threadIndex == 0) {
            shared.n_live_paths =
                (shared.n_live_paths < static_cast<int>(blockSize))
                    ? 0
                    : shared.n_live_paths - static_cast<int>(blockSize);
        }
        barrier.blockBarrier();
        if (has_path) {
            const unsigned int edge_pos = edge_size * path.x;
            const unsigned int nNei =
                d_output_graph[edge_pos + gbts_consts::nNei];
            const unsigned char level = d_levels[path.x];
            for (unsigned int nei = 0; nei < nNei; ++nei) {
                const unsigned int edge_idx =
                    d_output_graph[edge_pos + gbts_consts::nei_start + nei];
                if (level != d_levels[edge_idx] + 1) {
                    continue;
                }
                path_idx = vecmem::device_atomic_ref<unsigned int>(
                               nPathStoreSizeCounter)
                               .fetch_add(1u);
                if (path_idx >= nReachablePaths) {
                    break;
                }
                const int live_idx =
                    vecmem::device_atomic_ref<
                        int, vecmem::device_address_space::local>(
                        shared.n_live_paths)
                        .fetch_add(1);
                if (live_idx >=
                    static_cast<int>(
                        traccc::device::gbts_consts::live_path_buffer)) {
                    break;
                }
                d_path_store[path_idx] =
                    int2{static_cast<int>(edge_idx), static_cast<int>(path.y)};
                shared_live_paths[static_cast<unsigned int>(live_idx)] =
                    uint2{edge_idx, path_idx};
            }
        }
        barrier.blockBarrier();
    }
}

}  // namespace traccc::device
