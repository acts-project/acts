// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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

/// (Global Event Data) Payload for the @c traccc::device::gbts_fill_seeds_from_sorted_bids
/// function
struct gbts_fill_seeds_from_sorted_bids_payload {
  unsigned int nSeeds;
  /// seed quality (30 bits) + dropout_index (2 bits) + prop_idx (32 bits) for
  /// ordering the seeds
  vecmem::data::vector_view<unsigned long long int> seed_bids;
  /// Per-path (edge index, parent path-store index or -1) entries
  vecmem::data::vector_view<const int2> path_store;
  /// (node1, node2) original spacepoint indices per compacted edge
  vecmem::data::vector_view<const uint2> output_edge_nodes;
  /// Output: 3-SP seeds placed in this non-resizable buffer
  edm::seed_collection::view output_seeds;
};

/// @brief From an array of sorted bids convert accepted seed proposals into 3-spacepoint edm::seeds.
///
/// Each thread (strided by gridSize) processes one accepted proposal,
/// reads the three constituent spacepoint and places the seed in the
/// output buffer
///
/// @param[in] thread_id Thread identifier for the kernel launch
/// @param[in] payload   The global memory payload
///
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fill_seeds_from_sorted_bids(
    const thread_id_t& thread_id,
    const gbts_fill_seeds_from_sorted_bids_payload& payload);

}  // namespace traccc::device

#include "traccc/gbts_seeding/device/impl/gbts_fill_seeds_from_sorted_bids.ipp"
