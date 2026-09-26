/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/gbts_seeding/device/details/gbts_graph_edges_helpers.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_count_graph_edges(
    const thread_id_t& thread_id,
    const gbts_count_graph_edges_payload& payload) {
  const vecmem::device_vector<const float4> d_node_params(payload.node_params);
  const vecmem::device_vector<const float> d_node_phi(payload.node_phi);
  const vecmem::device_vector<const unsigned int> d_pair_work_begin(
      payload.pair_work_begin);
  const vecmem::device_vector<const gbts_graph_edges_work_item> d_work_items(
      payload.work_items);
  vecmem::device_vector<unsigned int> d_edge_counts(payload.edge_counts);
  vecmem::device_vector<unsigned int> d_num_outgoing_edges(
      payload.num_outgoing_edges);

  constexpr unsigned int chunk_size = gbts_consts::edge_chunk_size;
  const unsigned int threadIndex = thread_id.getLocalThreadIdX();
  const unsigned int nWork = d_pair_work_begin.back();

  for (unsigned int work = thread_id.getBlockIdX(); work < nWork;
       work += thread_id.getGridDimX()) {
    const gbts_graph_edges_work_item item = d_work_items[work];
    unsigned int count = 0u;
    if (threadIndex < item.chunk_size) {
      const unsigned int node1 = item.chunk_begin + threadIndex;
      // The lambda is called by the window-walking function for each edge
      // to count. It does not use the node params, but they are passed
      // for consistency with the fill kernel.
      auto count_edge = [&](const unsigned int, const float4&, const float,
                            const float4&) { ++count; };
      details::gbts_walk_window(
          d_node_phi, d_node_params, item.outer_begin, item.outer_end,
          d_node_params[node1], d_node_phi[node1], item.delta_phi,
          payload.gbts_make_graph_edges_params, count_edge);
      if (count > 0u) {
        vecmem::device_atomic_ref<unsigned int>(
            d_num_outgoing_edges[node1 + 1u])
            .fetch_add(count);
      }
    }
    d_edge_counts[work * chunk_size + threadIndex] = count;
  }
}

}  // namespace traccc::device
