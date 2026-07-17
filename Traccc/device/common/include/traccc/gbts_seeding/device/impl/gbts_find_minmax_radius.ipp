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
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_find_minmax_radius(
    const thread_id_t& thread_id,
    const gbts_find_minmax_radius_payload& payload) {

    const vecmem::device_vector<const unsigned int> d_eta_bin_views(
        payload.eta_bin_views);
    const vecmem::device_vector<const float4> d_node_params(
        payload.node_params);
    vecmem::device_vector<float> d_bin_rads(payload.bin_rads);

    const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
    const unsigned int blockDimX = thread_id.getBlockDimX();
    const unsigned int gridDimX = thread_id.getGridDimX();

    for (unsigned int globalIndex = globalIdx; globalIndex < payload.nEtaBins;
         globalIndex += blockDimX * gridDimX) {

        const unsigned int node_start = d_eta_bin_views[2u * globalIndex];
        const unsigned int node_end = d_eta_bin_views[2u * globalIndex + 1u];

        float min_r = 1e8f;
        float max_r = -1e8f;

        for (unsigned int node_idx = node_start; node_idx < node_end;
             node_idx++) {
            const float r = d_node_params[node_idx].z;
            max_r = fmaxf(r, max_r);
            min_r = fminf(r, min_r);
        }

        d_bin_rads[2u * globalIndex] = min_r;
        d_bin_rads[2u * globalIndex + 1u] = max_r;
    }
}

}  // namespace traccc::device
