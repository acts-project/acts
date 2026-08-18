/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device::details {

/// @brief Register a new seed proposal and let it bid for every edge on its
/// path.
///
/// Records the (quality, path index) tuple at prop_idx, packs the seed
/// bid (quality in high 32 bits, proposal id in low 32 bits), and walks the
/// path-store chain atomically claiming each edge it touches.  Whenever the
/// new bid loses, the proposal's own ambiguity tag is set to -1; whenever it
/// wins, the previous holder's ambiguity tag is set to -1.
///
/// @param[in]  qual                       Fit quality (higher wins)
/// @param[in]  path_idx                   Leaf entry in d_path_store_view
/// @param[in]  prop_idx                   Slot for the new proposal
/// @param[in,out] d_seed_ambiguity_view   Per-seed ambiguity tag
/// @param[in,out] d_seed_proposals_view   Per-seed (quality, path index)
/// @param[in,out] d_edge_bids_view        Per-edge highest-bidder seed
/// @param[in]  d_path_store_view          Per-path (parent, edge) entries
/// @param[in]  depth                      Optional maximum walk depth (-1 ==
/// unlimited)
///
TRACCC_HOST_DEVICE
inline void gbts_create_seed_candidate(
    const int qual, const int path_idx, const unsigned int prop_idx,
    const vecmem::data::vector_view<char> d_seed_ambiguity_view,
    const vecmem::data::vector_view<int2> d_seed_proposals_view,
    const vecmem::data::vector_view<unsigned long long int> d_edge_bids_view,
    const vecmem::data::vector_view<const int2>& d_path_store_view,
    char depth = -1);

}  // namespace traccc::device::details

#include "traccc/gbts_seeding/device/impl/gbts_create_seed_candidate.ipp"
