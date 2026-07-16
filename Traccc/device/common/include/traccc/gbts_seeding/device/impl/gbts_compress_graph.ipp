/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_compress_graph(
    const thread_id_t& thread_id, const gbts_compress_graph_payload& payload) {

    const vecmem::device_vector<const unsigned int> d_orig_node_index(
        payload.orig_node_index);
    const vecmem::device_vector<const uint2> d_edge_nodes(payload.edge_nodes);
    const vecmem::device_vector<const unsigned char> d_num_neighbours(
        payload.num_neighbours);
    const vecmem::device_vector<const unsigned int> d_neighbours(
        payload.neighbours);
    const vecmem::device_vector<const int> d_reIndexer(payload.reIndexer);
    vecmem::device_vector<unsigned int> d_output_graph(payload.output_graph);

    const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
    const unsigned int blockDimX = thread_id.getBlockDimX();
    const unsigned int gridDimX = thread_id.getGridDimX();

    for (unsigned int globalIndex = globalIdx; globalIndex < payload.nEdges;
         globalIndex += blockDimX * gridDimX) {

        const int newIdx = d_reIndexer[globalIndex];
        if (newIdx == -1) {
            continue;
        }

        // Row-major output graph: each edge owns a contiguous block of
        // edge_size = 2 + 1 + nMaxNei ints ([node1, node2, nNei,
        // nei0..neiN-1]).
        const unsigned int edge_size = 2u + 1u + payload.nMaxNei;
        const unsigned int pos = edge_size * static_cast<unsigned int>(newIdx);

        const uint2 edge_nodes = d_edge_nodes[globalIndex];
        d_output_graph[pos + gbts_consts::node1] =
            d_orig_node_index[edge_nodes.x];
        d_output_graph[pos + gbts_consts::node2] =
            d_orig_node_index[edge_nodes.y];

        const unsigned char nNei = d_num_neighbours[globalIndex];
        d_output_graph[pos + gbts_consts::nNei] = nNei;
        const unsigned int nei_pos = payload.nMaxNei * globalIndex;
        for (unsigned int k = 0u; k < nNei; k++) {
            d_output_graph[pos + gbts_consts::nei_start + k] =
                static_cast<unsigned int>(
                    d_reIndexer[d_neighbours[nei_pos + k]]);
        }
    }
}

}  // namespace traccc::device
