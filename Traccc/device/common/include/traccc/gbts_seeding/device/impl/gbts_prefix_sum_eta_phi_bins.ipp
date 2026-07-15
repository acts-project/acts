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

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_prefix_sum_eta_phi_bins(
    const thread_id_t& thread_id,
    const gbts_prefix_sum_eta_phi_bins_payload& payload) {

    const vecmem::device_vector<const unsigned int> d_eta_node_counter(
        payload.eta_node_counter);
    vecmem::device_vector<unsigned int> d_phi_cusums(payload.phi_cusums);

    const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
    const unsigned int blockDimX = thread_id.getBlockDimX();
    const unsigned int gridDimX = thread_id.getGridDimX();

    for (unsigned int globalIndex = globalIdx; globalIndex < payload.nEtaBins;
         globalIndex += blockDimX * gridDimX) {

        if (globalIndex == 0) {
            continue;
        }

        const unsigned int offset = payload.nPhiBins * globalIndex;
        const unsigned int val0 = d_eta_node_counter[globalIndex - 1];

        for (unsigned int phiIdx = 0; phiIdx < payload.nPhiBins; phiIdx++) {
            d_phi_cusums[offset + phiIdx] += val0;
        }
    }
}

}  // namespace traccc::device
