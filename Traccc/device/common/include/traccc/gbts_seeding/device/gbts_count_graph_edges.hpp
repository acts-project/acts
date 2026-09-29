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

/// Largest launch grid of gbts_count_graph_edges (the blocks stride over the
/// work items)
inline constexpr unsigned int gbts_count_graph_edges_max_blocks = 4096u;

/// (Global Event Data) Payload for the @c traccc::device::gbts_count_graph_edges
/// function
struct gbts_count_graph_edges_payload {
  /// Upper bound of the number of work items (sizes the launch grid)
  unsigned int nWorkMax;
  /// Per bin pair, first work item of the pair.
  vecmem::data::vector_view<const unsigned int> pair_work_begin;
  /// The work items, in bin pair order (from gbts_build_edge_work_list)
  vecmem::data::vector_view<const gbts_graph_edges_work_item> work_items;
  /// Per node, (tau_min, tau_max, r, z)
  vecmem::data::vector_view<const float4> node_params;
  /// Per node, phi
  vecmem::data::vector_view<const float> node_phi;
  /// Edge cuts
  traccc::gbts_make_graph_edges_params gbts_make_graph_edges_params;
  /// Output: edges per (work item, thread), at [work * chunk_size + thread]
  vecmem::data::vector_view<unsigned int> edge_counts;
  /// Output: edges per inner node, accumulated at [node + 1].
  vecmem::data::vector_view<unsigned int> num_outgoing_edges;
};

/// @brief Count the candidate edges of every inner node.
///
/// A work item is a bin pair and one chunk of the pair's inner bin. Every
/// thread owns one inner node of the chunk and walks its delta-phi window in
/// the phi-sorted outer bin, testing the candidates against the edge cuts.
///
/// @param[in] thread_id      Thread identifier (one block per work item)
/// @param[in,out] payload    The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_count_graph_edges(
    const thread_id_t& thread_id,
    const gbts_count_graph_edges_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_count_graph_edges.ipp"
