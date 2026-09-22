/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
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
TRACCC_HOST_DEVICE inline void gbts_finish_cca(
    const thread_id_t& thread_id, const gbts_finish_cca_payload& payload) {
  const vecmem::device_vector<const unsigned int> d_output_graph(
      payload.output_graph);
  const vecmem::device_vector<const unsigned char> d_levels(payload.levels);
  vecmem::device_vector<int2> d_outgoing_paths(payload.outgoing_paths);
  vecmem::device_vector<unsigned char> d_has_parent(payload.has_parent);

  // Row-major output graph: each edge owns [node1, node2, nNei, nei0..].
  const unsigned int edge_size =
      gbts_consts::nei_start + payload.max_num_neighbours;

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  for (unsigned int globalIndex = globalIdx;
       globalIndex < payload.nConnectedEdges;
       globalIndex += blockDimX * gridDimX) {
    const unsigned char level = d_levels[globalIndex];
    if (level > gbts_consts::max_seed_candidate_length) {
      // Longest path exceeds max_seed_candidate_length edges.
      d_outgoing_paths[globalIndex] = int2{0, -1};
      continue;
    }
    const unsigned int edge_pos = edge_size * globalIndex;
    const unsigned int nNei = d_output_graph[edge_pos + gbts_consts::nNei];
    for (unsigned int k = 0u; k < nNei; ++k) {
      d_has_parent[d_output_graph[edge_pos + gbts_consts::nei_start + k]] = 1u;
    }
    d_outgoing_paths[globalIndex].y =
        static_cast<int>(level >= payload.minLevel) - 1;
  }
}

}  // namespace traccc::device
