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

/// (Global Event Data) Payload for the @c traccc::device::gbts_reset_edge_bids
/// function
struct gbts_reset_edge_bids_payload {
  /// Number of seed proposals
  unsigned int nProps;
  /// Per-path (edge index, parent path-store index or -1) entries
  vecmem::data::vector_view<const int2> path_store;
  /// In/out: per-seed-proposal (path_store index, level)
  vecmem::data::vector_view<int2> seed_proposals;
  /// In/out: per-edge highest-bidder seed proposal (cleared between
  /// rounds)
  vecmem::data::vector_view<unsigned long long int> edge_bids;
  /// In/out: per-seed-proposal ambiguity tag
  vecmem::data::vector_view<char> seed_ambiguity;
  /// In/out: global atomic counter of rejected proposals
  unsigned int* nRejectedPropsCounter;
};

/// @brief Mark a losing seed proposal against the current edge bids.
///
/// Processes one proposal (the grid-stride loop lives in the kernel wrapper):
/// compares it against the winning bid recorded in the edge bids, and either
/// updates the ambiguity tag or atomically increments the rejected counter
/// if it loses.
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_reset_edge_bids(
    const thread_id_t& thread_id, const gbts_reset_edge_bids_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_reset_edge_bids.ipp"
