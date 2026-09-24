/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
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
  const vecmem::device_vector<const unsigned int> d_num_outgoing_edges(
      payload.num_outgoing_edges);
  const vecmem::device_vector<const unsigned int> d_reIndexer(
      payload.reIndexer);
  vecmem::device_vector<unsigned int> d_output_graph(payload.output_graph);

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  const unsigned int nEdgesTotal = d_num_outgoing_edges.back();
  const unsigned int nEdges = math::min(nEdgesTotal, payload.nEdgesMax);
  for (unsigned int globalIndex = globalIdx; globalIndex < nEdges;
       globalIndex += blockDimX * gridDimX) {
    const unsigned int scan = d_reIndexer[globalIndex];
    const unsigned int prev =
        (globalIndex == 0u) ? 0u : d_reIndexer[globalIndex - 1u];
    if (scan == prev) {
      continue;
    }
    const unsigned int newIdx = scan - 1u;
    if (newIdx >= payload.nConnectedEdgesMax) {
      continue;
    }

    // Row-major output graph: each edge owns a contiguous block of
    // nei_start + nMaxNei ints ([node1, node2, nNei, nei0..neiN-1]).
    const unsigned int edge_size = gbts_consts::nei_start + payload.nMaxNei;
    const unsigned int pos = edge_size * newIdx;

    const uint2 edge_nodes = d_edge_nodes[globalIndex];
    d_output_graph[pos + gbts_consts::node1] = d_orig_node_index[edge_nodes.x];
    d_output_graph[pos + gbts_consts::node2] = d_orig_node_index[edge_nodes.y];

    const unsigned char nNei = d_num_neighbours[globalIndex];
    const unsigned int nei_pos = payload.nMaxNei * globalIndex;
    // Neighbours beyond the capacity were dropped above and are dropped
    // from the row, so every index in the compacted graph is in bounds.
    unsigned int kept = 0u;
    for (unsigned int k = 0u; k < nNei; k++) {
      // Every recorded neighbour is itself kept.
      const unsigned int nei = d_reIndexer[d_neighbours[nei_pos + k]] - 1u;
      if (nei >= payload.nConnectedEdgesMax) {
        continue;
      }
      d_output_graph[pos + gbts_consts::nei_start + kept] = nei;
      ++kept;
    }
    d_output_graph[pos + gbts_consts::nNei] = kept;
  }
}

}  // namespace traccc::device
