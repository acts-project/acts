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

// System include(s).
#include <cstdint>

namespace traccc::device {

/// (Global Event Data) Payload for the @c
/// traccc::device::gbts_bid_seeds_for_hits function
struct gbts_bid_seeds_for_hits_payload {
    /// Number of seed proposals
    unsigned int nProps;
    /// Number of accepted seeds (nProps - nRejectedProps)
    unsigned int nSeeds;
    /// Per-edge row stride in the output graph (= 2 + 1 +
    /// max_num_neighbours)
    unsigned int edge_size;
    /// Compacted graph from gbts_compress_graph
    vecmem::data::vector_view<const unsigned int> output_graph;
    /// Per-seed-proposal (path_store index, level)
    vecmem::data::vector_view<const int2> seed_proposals;
    /// Per-path (edge index, parent path-store index or -1) entries
    vecmem::data::vector_view<const int2> path_store;
    /// Per-seed-proposal ambiguity tag
    vecmem::data::vector_view<const char> seed_ambiguity;
    /// In/out: per-hit highest-bidder seed (packed 64-bit)
    vecmem::data::vector_view<unsigned long long int> hit_bids;
};

/// @brief One accepted seed bids on its constituent hits.
///
/// Processes one proposal (the grid-stride loop lives in the kernel wrapper):
/// walks the seed's edges via the compact graph to enumerate hit indices, and
/// atomically updates the hit bids with its packed seed bid if it outranks
/// the current best for that hit.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_bid_seeds_for_hits(
    const thread_id_t& thread_id,
    const gbts_bid_seeds_for_hits_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_bid_seeds_for_hits.ipp"
