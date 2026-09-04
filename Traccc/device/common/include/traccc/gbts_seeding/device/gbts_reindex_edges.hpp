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

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_reindex_edges
/// function
struct gbts_reindex_edges_payload {
  /// Number of original edges
  unsigned int nEdges;
  /// In/out: per-edge "kept" flag in, compact new index out
  vecmem::data::vector_view<int> reIndexer;
  /// In/out: global atomic counter of edges that survived re-indexing
  unsigned int* nConnectedEdgesCounter;
};

/// @brief Replace the per-edge "kept" flag with its compacted index.
///
/// Each thread reads its slot in the re-indexer; if the edge is marked
/// "kept", it atomically claims the next compact slot and writes that
/// compact index back; otherwise the slot is set to a sentinel.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_reindex_edges(
    const thread_id_t& thread_id, const gbts_reindex_edges_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_reindex_edges.ipp"
