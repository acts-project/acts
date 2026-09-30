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
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_finish_cca
/// function
struct gbts_finish_cca_payload {
  /// Number of edges in the compacted graph
  unsigned int nConnectedEdges;
  /// Maximum number of neighbours retained per edge
  unsigned int max_num_neighbours;
  /// Minimum level (path length in edges) of a seed root
  unsigned char minLevel;
  /// Compacted graph from gbts_compress_graph
  vecmem::data::vector_view<const unsigned int> output_graph;
  /// Per-edge level from gbts_run_cca_iteration
  vecmem::data::vector_view<const unsigned char> levels;
  /// In/out: per-edge (subtree path count, terminus flag: 0 = seed root
  /// candidate, -1 = too short or not settled)
  vecmem::data::vector_view<int2> outgoing_paths;
  /// Output: per-edge "has a parent" mark, zero-initialised
  vecmem::data::vector_view<unsigned char> has_parent;
};

/// @brief Write the parent marks and terminus flags of the settled edges.
///
/// Run once after the sweeps of gbts_run_cca_iteration. Every settled edge
/// marks its neighbours as having a parent and gets its terminus flag from
/// its level. An edge whose longest path exceeds max_seed_candidate_length
/// edges is not a root and marks nobody, which truncates the path.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_finish_cca(
    const thread_id_t& thread_id, const gbts_finish_cca_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_finish_cca.ipp"
