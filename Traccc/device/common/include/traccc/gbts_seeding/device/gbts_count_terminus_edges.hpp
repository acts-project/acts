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
/// traccc::device::gbts_count_terminus_edges function
struct gbts_count_terminus_edges_payload {
    /// Number of edges in the compacted graph
    unsigned int nConnectedEdges;
    /// Per-edge longest-outgoing-path summary from CCA
    vecmem::data::vector_view<short2> outgoing_paths;
    /// Total number of paths reachable from any terminus edge
    unsigned int* nPathsCounter;
    /// Running size of the path store (initialised to nTerminusEdges)
    unsigned int* nPathStoreSizeCounter;
};

/// @brief Count terminus edges (those with no live outgoing path) and total
/// paths.
///
/// Each thread inspects one edge; terminus edges atomically claim a slot in
/// the path store and fold their reachable-path count into the total.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_count_terminus_edges(
    const thread_id_t& thread_id,
    const gbts_count_terminus_edges_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_count_terminus_edges.ipp"
