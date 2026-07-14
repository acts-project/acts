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
#include "traccc/edm/seed_collection.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Global Event Data) Payload for the @c traccc::device::gbts_convert_seeds
/// function
struct gbts_convert_seeds_payload {
    /// Number of seed proposals
    unsigned int nProps;
    /// Number of accepted seeds (nProps - nRejectedProps)
    unsigned int nSeeds;
    /// Maximum number of neighbours retained per edge
    unsigned int max_num_neighbours;
    /// Per-seed-proposal (path_store index, level)
    vecmem::data::vector_view<const int2> seed_proposals;
    /// Per-seed-proposal ambiguity tag
    vecmem::data::vector_view<const char> seed_ambiguity;
    /// Per-path (edge index, parent path-store index or -1) entries
    vecmem::data::vector_view<const int2> path_store;
    /// Compacted graph from gbts_compress_graph
    vecmem::data::vector_view<const unsigned int> output_graph;
    /// Reduced (x, y, z, r) per original spacepoint
    vecmem::data::vector_view<const float4> reducedSP;
    /// Output: 3-SP seeds appended to this resizable buffer
    edm::seed_collection::view output_seeds;
    /// Per-hit highest-bidder seed (read for dropout decisions)
    vecmem::data::vector_view<unsigned long long int> hit_bids;
    /// Dropout / curvature / ambiguity cut parameters
    traccc::gbts_convert_seeds_params gbts_convert_seeds_params;
};

/// @brief Convert accepted seed proposals into 3-spacepoint edm::seeds.
///
/// Each thread (strided by gridSize) processes one accepted proposal,
/// reads the three constituent spacepoints, applies dropout / curvature /
/// hit-fraction cuts, and appends a seed to the output resizable buffer on
/// success.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_convert_seeds(
    const thread_id_t& thread_id, const gbts_convert_seeds_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_convert_seeds.ipp"
