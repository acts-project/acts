// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/gbts_seeding/device/details/gbts_seed_conversion_helpers.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <array>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fill_seeds_from_sorted_bids(
    const thread_id_t& thread_id,
    const gbts_fill_seeds_from_sorted_bids_payload& payload) {
  edm::seed_collection::device seeds_device(payload.output_seeds);
  vecmem::device_vector<long long unsigned int> d_seed_bids(payload.seed_bids);
  const vecmem::device_vector<const int2> d_path_store(payload.path_store);
  const vecmem::device_vector<const uint2> d_output_edge_nodes(
      payload.output_edge_nodes);

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  for (unsigned int prop_idx = globalIdx; prop_idx < payload.nSeeds;
       prop_idx += blockDimX * gridDimX) {
    // take the dropout index from the saved bid
    unsigned long long int bid = d_seed_bids.at(prop_idx);
    int2 prop = {static_cast<int>(bid >> 32),
                 static_cast<int>(bid & 0xFFFFFFFF)};
    unsigned int dropout_index = static_cast<unsigned int>(prop.x & 0x3);
    // match quality to non-sorted case
    prop.x -= dropout_index;

    details::Tracklet seed;
    seed.size = 0;
    int2 path = int2{0, prop.y};
    while (path.y >= 0) {
      path = d_path_store[static_cast<unsigned int>(path.y)];
      seed.nodes[seed.size++] =
          d_output_edge_nodes[static_cast<unsigned int>(path.x)].x;
    }
    seed.nodes[seed.size++] =
        d_output_edge_nodes[static_cast<unsigned int>(path.x)].y;

    const int node1 = seed.nodes[seed.size - 1 - (dropout_index == 0)];
    const int node2 = seed.nodes[(seed.size - 1) / 2 + (dropout_index == 2)];
    const int node3 = seed.nodes[0];

    auto out = seeds_device.at(prop_idx);
    details::set_seed(out, node1, node2, node3, static_cast<float>(prop.x));
  }
}

}  // namespace traccc::device
