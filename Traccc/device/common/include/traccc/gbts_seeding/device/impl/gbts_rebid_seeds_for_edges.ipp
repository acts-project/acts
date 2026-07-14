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
#include "traccc/gbts_seeding/device/details/gbts_create_seed_candidate.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_rebid_seeds_for_edges(
    const thread_id_t& thread_id,
    const gbts_rebid_seeds_for_edges_payload& payload) {

    vecmem::device_vector<char> d_seed_ambiguity(payload.seed_ambiguity);
    vecmem::device_vector<int2> d_seed_proposals(payload.seed_proposals);

    const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
    const unsigned int blockDimX = thread_id.getBlockDimX();
    const unsigned int gridDimX = thread_id.getGridDimX();

    for (unsigned int prop_idx = globalIdx; prop_idx < payload.nProps;
         prop_idx += blockDimX * gridDimX) {

        const char ambi = d_seed_ambiguity[prop_idx];

        if (payload.first_round) {
            if (ambi == 0) {
                // rebid 'best seed from edge' in later rounds
                d_seed_ambiguity[prop_idx] = 1;
                // Here there is no return by design
            } else {
                d_seed_ambiguity[prop_idx] = -2;
                // count rejected props to calculate nSeeds
                vecmem::device_atomic_ref<unsigned int>(
                    *payload.nRejectedPropsCounter)
                    .fetch_add(1u);
                continue;
            }
        } else if ((ambi == -2) | (ambi == 0)) {
            // only rebid for maybes
            continue;
        }
        const int2 prop = d_seed_proposals[prop_idx];

        details::gbts_create_seed_candidate(
            prop.x, prop.y, prop_idx, payload.seed_ambiguity,
            payload.seed_proposals, payload.edge_bids, payload.path_store, -1);
    }
}

}  // namespace traccc::device
