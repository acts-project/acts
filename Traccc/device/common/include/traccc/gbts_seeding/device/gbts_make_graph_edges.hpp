/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
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

/// (Global Event Data) Payload for the @c traccc::device::gbts_make_graph_edges
/// function
struct gbts_make_graph_edges_payload {
    /// Number of bin-pair tasks (also the CUDA block count)
    unsigned int nUsedBinPairs;
    /// Upper bound on the number of edges to write
    unsigned int nMaxEdges;
    /// Number of phi bins per eta slice
    unsigned int nPhiBins;
    /// Per-bin-pair (begin1, end1, begin2, end2) node ranges, flat
    vecmem::data::vector_view<const unsigned int> bin_pair_views;
    /// Per-bin-pair max delta-phi window for edge candidates
    vecmem::data::vector_view<const float> bin_pair_dphi;
    /// Per-node (tau_min, tau_max, r, z)
    vecmem::data::vector_view<const float4> node_params;
    /// Per-node phi
    vecmem::data::vector_view<const float> node_phi;
    /// Edge-making geometric / kinematic cuts
    traccc::gbts_make_graph_edges_params gbts_make_graph_edges_params;
    /// In/out: global atomic counter for the next edge slot to write
    unsigned int* nEdgesCounter;
    /// Output: (src, dst) node indices per edge
    vecmem::data::vector_view<uint2> edge_nodes;
    /// Output: packed per-edge [exp(-eta), curv, phi_z, phi_w] used by
    /// matching
    vecmem::data::vector_view<float4> edge_params;
    /// Output: per-destination-node incoming-edge count (atomic)
    vecmem::data::vector_view<unsigned int> num_outgoing_edges;
};

/// (Shared Event Data) Payload for the @c traccc::device::gbts_make_graph_edges
/// function
///
/// Shared-memory scratch for gbts_make_graph_edges: a block-local copy of the
/// current node1 chunk (phi values and packed node params).
struct gbts_make_graph_edges_shared_payload {
    /// Shared-mem cache: phi / node
    vecmem::data::vector_view<float> phi;
    /// Shared-mem cache: (tau_min, tau_max, r, z) float4 / node
    vecmem::data::vector_view<float4> node_pack;
};

/// @brief Create candidate edges between node pairs in compatible (eta, phi)
/// bins.
///
/// One CUDA block handles one bin-pair task.  Threads stage a chunk of bin-1
/// nodes into shared-memory caches (phi as a separate float array, and the
/// (tau_min, tau_max, r, z) float4 per node), then for every bin-2 node test
/// the cached chunk against geometric and kinematic cuts
/// (gbts_check_edge_candidate), atomically reserving a slot in the output via
/// the edge counter.
///
/// @param[in] thread_id          Thread/block identifier (one block/task)
/// @param[in] barrier            Block-wide barrier
/// @param[in,out] payload        The global memory payload
/// @param[in,out] shared_payload The shared memory payload
///
template <concepts::thread_id1 thread_id_t, concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void gbts_make_graph_edges(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const gbts_make_graph_edges_payload& payload,
    const gbts_make_graph_edges_shared_payload& shared_payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_make_graph_edges.ipp"
