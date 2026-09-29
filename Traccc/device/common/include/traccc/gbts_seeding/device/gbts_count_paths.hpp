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
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Global Event Data) Payload for the @c
/// traccc::device::gbts_count_paths function
struct gbts_count_paths_payload {
  /// Number of edges in the compacted graph
  unsigned int nConnectedEdges;
  /// Per-edge (subtree path count, terminus flag) from CCA
  vecmem::data::vector_view<const int2> outgoing_paths;
  /// Per-edge "has a settled parent" mark from CCA (roots have none)
  vecmem::data::vector_view<const unsigned char> has_parent;
  /// Output: per-edge number of paths, 1 + subtree for a terminus edge and 0
  /// otherwise. The launcher scans it in place.
  vecmem::data::vector_view<unsigned int> path_counts;
};

/// @brief Count the paths starting at every terminus edge.
///
/// The launcher's prefix sum over the counts gives every terminus edge a
/// contiguous range of the path store, which gbts_fill_path_store fills.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_count_paths(
    const thread_id_t& thread_id, const gbts_count_paths_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_count_paths.ipp"
