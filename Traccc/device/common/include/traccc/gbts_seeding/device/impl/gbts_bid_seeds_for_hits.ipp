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
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_bid_seeds_for_hits(
    const thread_id_t& thread_id,
    const gbts_bid_seeds_for_hits_payload& payload) {

    const vecmem::device_vector<const unsigned int> d_output_graph(
        payload.output_graph);
    const vecmem::device_vector<const int2> d_path_store(payload.path_store);
    const vecmem::device_vector<const char> d_seed_ambiguity(
        payload.seed_ambiguity);
    const vecmem::device_vector<const int2> d_seed_proposals(
        payload.seed_proposals);
    vecmem::device_vector<unsigned long long int> d_hit_bids(payload.hit_bids);

    const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
    const unsigned int blockDimX = thread_id.getBlockDimX();
    const unsigned int gridDimX = thread_id.getGridDimX();

    for (unsigned int prop_idx = globalIdx; prop_idx < payload.nProps;
         prop_idx += blockDimX * gridDimX) {

        if (d_seed_ambiguity[prop_idx] == -2) {
            continue;
        }
        const int2 prop = d_seed_proposals[prop_idx];
        const unsigned long long int seed_bid =
            (static_cast<unsigned long long int>(prop.x) << 32) |
            (static_cast<unsigned long long int>(prop_idx));

        int2 path = int2{0, prop.y};
        while (path.y >= 0) {
            path = d_path_store[static_cast<unsigned int>(path.y)];
            const unsigned int sp_idx =
                d_output_graph[payload.edge_size *
                                   static_cast<unsigned int>(path.x) +
                               gbts_consts::node1];
            vecmem::device_atomic_ref<unsigned long long int> atomic_bid(
                d_hit_bids[sp_idx]);
            atomic_bid.fetch_max(seed_bid);
        }
        const unsigned int sp_idx =
            d_output_graph[payload.edge_size *
                               static_cast<unsigned int>(path.x) +
                           gbts_consts::node2];
        vecmem::device_atomic_ref<unsigned long long int> atomic_bid(
            d_hit_bids[sp_idx]);
        atomic_bid.fetch_max(seed_bid);
    }
}

}  // namespace traccc::device
