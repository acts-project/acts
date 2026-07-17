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
/// traccc::device::gbts_add_terminus_to_path_store function
struct gbts_add_terminus_to_path_store_payload {
    /// Number of edges in the compacted graph
    unsigned int nConnectedEdges;
    /// Output: per-path (edge index, parent path-store index or -1)
    /// entries; terminus rows occupy the first nTerminusEdges slots
    vecmem::data::vector_view<int2> path_store;
    /// Per-edge longest-outgoing-path summary from CCA
    vecmem::data::vector_view<const short2> outgoing_paths;
};

/// @brief Seed the path store with one entry per terminus edge.
///
/// Each thread inspects one edge; if it is a terminus (no live outgoing path),
/// it writes a (parent = -1, edge) record into the path store, marking the
/// root of a future path traversal.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_add_terminus_to_path_store(
    const thread_id_t& thread_id,
    const gbts_add_terminus_to_path_store_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_add_terminus_to_path_store.ipp"
