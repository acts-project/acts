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
TRACCC_HOST_DEVICE inline void gbts_run_cca_iteration(
    const thread_id_t& thread_id,
    const gbts_run_cca_iteration_payload& payload) {
  const vecmem::device_vector<const unsigned int> d_output_graph(
      payload.output_graph);
  vecmem::device_vector<unsigned char> d_levels(payload.levels);
  vecmem::device_vector<int2> d_outgoing_paths(payload.outgoing_paths);
  vecmem::device_vector<unsigned int> d_changed(payload.changed);

  // Local copies: the writes to the char-typed levels may alias the payload,
  // which would force a reload of its fields in every iteration.
  const unsigned int nConnectedEdges = payload.nConnectedEdges;
  const unsigned int sweep = payload.iter;
  const bool firstSweep = (sweep == 0u);

  // Nothing changed in the previous sweep: the levels are final.
  if (!firstSweep && (d_changed[sweep - 1u] == 0u)) {
    return;
  }

  // Level of an edge whose longest path exceeds max_seed_candidate_length
  // edges.
  constexpr unsigned int unsettledLevel =
      gbts_consts::max_seed_candidate_length + 1u;
  // Row-major output graph: each edge owns [node1, node2, nNei, nei0..].
  const unsigned int edge_size =
      gbts_consts::nei_start + payload.max_num_neighbours;

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  for (unsigned int globalIndex = globalIdx; globalIndex < nConnectedEdges;
       globalIndex += blockDimX * gridDimX) {
    // Walk the edges from the back to the front because the longest paths are
    // expected to be found in the last edges.
    const unsigned int edgeIdx = nConnectedEdges - 1u - globalIndex;
    const unsigned int edge_pos = edge_size * edgeIdx;
    const unsigned int nNeighbours =
        d_output_graph[edge_pos + gbts_consts::nNei];

    // The level is one more than the highest neighbour level.
    unsigned int maxNeighbourLevel = 0u;
    for (unsigned int k = 0u; k < nNeighbours; ++k) {
      const unsigned int neighbourIdx =
          d_output_graph[edge_pos + gbts_consts::nei_start + k];
      const unsigned int neighbourLevel = d_levels[neighbourIdx];
      if (neighbourLevel > maxNeighbourLevel) {
        maxNeighbourLevel = neighbourLevel;
      }
    }
    unsigned int level = 1u + maxNeighbourLevel;
    if (level > unsettledLevel) {
      level = unsettledLevel;
    }

    // Paths below the edge. One for every neighbour on a longest path, plus
    // the paths below that neighbour.
    int nPathsBelow = 0;
    if (level < unsettledLevel) {
      for (unsigned int k = 0u; k < nNeighbours; ++k) {
        const unsigned int neighbourIdx =
            d_output_graph[edge_pos + gbts_consts::nei_start + k];
        if (d_levels[neighbourIdx] + 1u == level) {
          nPathsBelow +=
              1 + (firstSweep ? 0 : d_outgoing_paths[neighbourIdx].x);
        }
      }
    }

    if (firstSweep || (d_levels[edgeIdx] != level) ||
        (d_outgoing_paths[edgeIdx].x != nPathsBelow)) {
      d_levels[edgeIdx] = static_cast<unsigned char>(level);
      d_outgoing_paths[edgeIdx].x = nPathsBelow;
      // Mark that the next sweep must be run.
      d_changed[sweep] = 1u;
    }
  }
}

}  // namespace traccc::device
