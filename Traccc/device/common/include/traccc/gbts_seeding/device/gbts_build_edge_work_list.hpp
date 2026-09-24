/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
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
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// One graph-making work item: one gbts_consts::edge_chunk_size-sized chunk
/// of the inner bin of one bin pair, resolved by
/// @c traccc::device::gbts_build_edge_work_list
struct gbts_graph_edges_work_item {
  unsigned int pair;         ///< bin pair
  unsigned int chunk;        ///< chunk index inside the inner bin
  unsigned int chunk_begin;  ///< first inner node of the chunk
  unsigned int chunk_size;   ///< inner nodes in the chunk
  unsigned int outer_begin;  ///< outer bin node range
  unsigned int outer_end;
  float delta_phi;  ///< delta-phi window of the bin pair
};

/// Block size of the (single block) gbts_build_edge_work_list kernel; must
/// be a power of two.
inline constexpr unsigned int gbts_build_edge_work_list_block_size = 1024u;

/// (Global Event Data) Payload for the @c
/// traccc::device::gbts_build_edge_work_list function
struct gbts_build_edge_work_list_payload {
  /// Number of bin pairs
  unsigned int nBinPairs;
  /// Per bin pair: (bin1, bin2)
  vecmem::data::vector_view<const uint2> bin_pairs;
  /// Per eta bin node offsets.
  vecmem::data::vector_view<const unsigned int> eta_bin_offsets;
  /// Per eta bin (min r, max r), flat (from gbts_find_minmax_radius)
  vecmem::data::vector_view<const float> bin_rads;
  /// Delta-phi window of a bin pair from its radial separation
  traccc::gbts_dphi_window_params gbts_dphi_window_params;
  /// Output: per bin pair the first graph-making work item of the pair.
  vecmem::data::vector_view<unsigned int> pair_work_begin;
  /// Output: the resolved work items, in bin pair order; sized for the
  /// host-side upper bound of the work item count
  vecmem::data::vector_view<gbts_graph_edges_work_item> work_items;
};

/// (Shared Event Data) Payload for the @c
/// traccc::device::gbts_build_edge_work_list function
struct gbts_build_edge_work_list_shared_payload {
  /// Block scan scratch, gbts_build_edge_work_list_block_size entries
  vecmem::data::vector_view<unsigned int> scratch;
};

/// @brief Lay out the graph-making work items on the device.
///
/// A work item is one bin pair and one gbts_consts::edge_chunk_size-sized
/// chunk of the pair's inner bin pairs. The items carry their node ranges
/// and delta-phi window, so the count and fill passes do not touch the
/// bin tables.
///
/// @param[in] thread_id      Thread identifier (one block)
/// @param[in] barrier        Block-wide barrier
/// @param[in,out] payload    The global memory payload
/// @param[in] shared_payload The shared memory payload
///
template <concepts::thread_id1 thread_id_t, concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void gbts_build_edge_work_list(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const gbts_build_edge_work_list_payload& payload,
    const gbts_build_edge_work_list_shared_payload& shared_payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_build_edge_work_list.ipp"
