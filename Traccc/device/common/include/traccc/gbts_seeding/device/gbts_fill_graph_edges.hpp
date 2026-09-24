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
#include "traccc/gbts_seeding/device/gbts_build_edge_work_list.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Largest launch grid of gbts_fill_graph_edges (the blocks stride over the
/// work items)
inline constexpr unsigned int gbts_fill_graph_edges_max_blocks = 16384u;

/// (Global Event Data) Payload for the @c traccc::device::gbts_fill_graph_edges
/// function
struct gbts_fill_graph_edges_payload {
  /// Upper bound of the number of work items (sizes the launch grid)
  unsigned int nWorkMax;
  /// Per bin pair, first work item of the pair.
  vecmem::data::vector_view<const unsigned int> pair_work_begin;
  /// The work items, in bin pair order (from gbts_build_edge_work_list)
  vecmem::data::vector_view<const gbts_graph_edges_work_item> work_items;
  /// Per node, (tau_min, tau_max, r, z)
  vecmem::data::vector_view<const float4> node_params;
  /// Per node, phi.
  vecmem::data::vector_view<const float> node_phi;
  /// Edge cuts
  traccc::gbts_make_graph_edges_params gbts_make_graph_edges_params;
  /// Per bin pair, first pair with the same inner bin
  vecmem::data::vector_view<const unsigned int> pair_group_begin;
  /// Edges per (work item, thread), from the count pass
  vecmem::data::vector_view<const unsigned int> edge_counts;
  /// Per node, [node] and [node + 1] are the begin and end of the node's
  /// edge bucket.
  vecmem::data::vector_view<const unsigned int> num_outgoing_edges;
  /// Packs the edge parameters
  edge_params_converter edge_params_maker;
  /// Capacity of the edge buffers; edges beyond it are dropped
  unsigned int nEdgesMax;
  /// Output: (outer node, inner node) per edge
  vecmem::data::vector_view<uint2> edge_nodes;
  /// Output: packed (eta, curvature, phi at node2, phi at node1) per edge
  vecmem::data::vector_view<short4> edge_params;
};

/// @brief Write the candidate edges in canonical order.
///
/// Same work decomposition as gbts_count_graph_edges. The edges of an inner
/// node form one bucket, ordered by bin pair and by outer node inside
/// a pair, so every thread knows the slot of its edges
///
///
/// @param[in] thread_id      Thread identifier (one block per work item)
/// @param[in,out] payload    The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fill_graph_edges(
    const thread_id_t& thread_id, const gbts_fill_graph_edges_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_fill_graph_edges.ipp"
