/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/gbts_seeding/device/details/gbts_graph_edges_helpers.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

namespace details {

/// Write the edge node1 -> node2 at @c cursor: the node pair (outer, inner)
/// and the packed parameters (eta, curvature, phi extrapolated to node2, phi
/// extrapolated to node1) used by the matching.
TRACCC_HOST_DEVICE inline void gbts_write_edge(
    const unsigned int cursor, const unsigned int node1,
    const unsigned int node2, const float4 np1, const float4 np2,
    const float phi1, const float phi2, const float4 geo,
    const gbts_make_graph_edges_params& cuts,
    const edge_params_converter& converter,
    vecmem::device_vector<uint2>& edge_nodes,
    vecmem::device_vector<short4>& edge_params) {
  const float eta = -math::log(math::sqrt(1.0f + geo.x * geo.x) - geo.x);
  const bool long_edge = (cuts.long_edge_dz < math::fabs(geo.z)) ||
                         (cuts.long_edge_dr < math::fabs(geo.y));
  edge_nodes[cursor] = uint2{node2, node1};
  edge_params[cursor] = converter.make_edge_params(
      eta, geo.w, phi2 + geo.w * np2.z, phi1 + geo.w * np1.z, long_edge);
}

}  // namespace details

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fill_graph_edges(
    const thread_id_t& thread_id,
    const gbts_fill_graph_edges_payload& payload) {
  const vecmem::device_vector<const float4> d_node_params(payload.node_params);
  const vecmem::device_vector<const float> d_node_phi(payload.node_phi);
  const vecmem::device_vector<const unsigned int> d_pair_work_begin(
      payload.pair_work_begin);
  const vecmem::device_vector<const gbts_graph_edges_work_item> d_work_items(
      payload.work_items);
  const vecmem::device_vector<const unsigned int> d_pair_group_begin(
      payload.pair_group_begin);
  const vecmem::device_vector<const unsigned int> d_edge_counts(
      payload.edge_counts);
  const vecmem::device_vector<const unsigned int> d_num_outgoing_edges(
      payload.num_outgoing_edges);
  vecmem::device_vector<uint2> d_edge_nodes(payload.edge_nodes);
  vecmem::device_vector<short4> d_edge_params(payload.edge_params);

  constexpr unsigned int chunk_size = gbts_consts::edge_chunk_size;
  const gbts_make_graph_edges_params& cuts =
      payload.gbts_make_graph_edges_params;
  const unsigned int threadIndex = thread_id.getLocalThreadIdX();
  const unsigned int nWork = d_pair_work_begin.back();

  for (unsigned int work = thread_id.getBlockIdX(); work < nWork;
       work += thread_id.getGridDimX()) {
    const gbts_graph_edges_work_item item = d_work_items[work];
    if (threadIndex >= item.chunk_size) {
      continue;
    }
    const unsigned int node1 = item.chunk_begin + threadIndex;
    const float4 np1 = d_node_params[node1];
    const float phi1 = d_node_phi[node1];

    // The edges of this work item start after those of the preceding pairs
    // of the same inner bin (same chunk index, same thread).
    unsigned int cursor = d_num_outgoing_edges[node1];
    for (unsigned int p = d_pair_group_begin[item.pair]; p < item.pair; ++p) {
      const unsigned int first = d_pair_work_begin[p];
      if (d_pair_work_begin[p + 1u] > first) {
        cursor +=
            d_edge_counts[(first + item.chunk) * chunk_size + threadIndex];
      }
    }
    const unsigned int cursor_end =
        math::min(d_num_outgoing_edges[node1 + 1u], payload.nEdgesMax);
    auto write_edge = [&](const unsigned int node2, const float4 np2,
                          const float phi2, const float4 geo) {
      if (cursor < cursor_end) {
        details::gbts_write_edge(cursor, node1, node2, np1, np2, phi1, phi2,
                                 geo, cuts, payload.edge_params_maker,
                                 d_edge_nodes, d_edge_params);
        ++cursor;
      }
    };

    details::gbts_walk_window(d_node_phi, d_node_params, item.outer_begin,
                              item.outer_end, np1, phi1, item.delta_phi, cuts,
                              write_edge);
  }
}

}  // namespace traccc::device
