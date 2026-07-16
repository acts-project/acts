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

/// (Global Event Data) Payload for the @c traccc::device::gbts_link_graph_edges
/// function
struct gbts_link_graph_edges_payload {
    /// Number of edges produced earlier
    unsigned int nEdges;
    /// (src, dst) node indices per edge
    vecmem::data::vector_view<const uint2> edge_nodes;
    /// Output: per-edge slot in the per-node incoming-edge list
    vecmem::data::vector_view<unsigned int> edge_links;
    /// In/out: per-node prefix-sum / write cursor of incoming edges
    vecmem::data::vector_view<unsigned int> num_outgoing_edges;
};

/// @brief Compute each edge's slot in its destination node's incoming-edge
/// list.
///
/// One thread per edge atomically increments the per-destination-node count
/// and records the returned slot in the edge-link list.  After this kernel
/// the count buffer has been turned into a write cursor that
/// gbts_match_graph_edges can read sequentially.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_link_graph_edges(
    const thread_id_t& thread_id, const gbts_link_graph_edges_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_link_graph_edges.ipp"
