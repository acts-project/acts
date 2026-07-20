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
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_link_graph_edges(
    const thread_id_t& thread_id,
    const gbts_link_graph_edges_payload& payload) {
  const vecmem::device_vector<const uint2> d_edge_nodes(payload.edge_nodes);
  vecmem::device_vector<unsigned int> d_edge_links(payload.edge_links);
  vecmem::device_vector<unsigned int> d_num_outgoing_edges(
      payload.num_outgoing_edges);

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  for (unsigned int globalIndex = globalIdx; globalIndex < payload.nEdges;
       globalIndex += blockDimX * gridDimX) {
    const unsigned int sharedNode = d_edge_nodes[globalIndex].y;
    const unsigned int pos = vecmem::device_atomic_ref<unsigned int>(
                                 d_num_outgoing_edges[sharedNode])
                                 .fetch_sub(1u);
    d_edge_links[pos - 1u] = static_cast<unsigned int>(globalIndex);
  }
}

}  // namespace traccc::device
