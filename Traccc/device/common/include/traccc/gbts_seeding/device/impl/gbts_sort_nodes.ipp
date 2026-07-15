/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_sort_nodes(
    const thread_id_t& thread_id, const gbts_sort_nodes_payload& payload) {

    const vecmem::device_vector<const float4> d_sp_params(payload.sp_params);
    const vecmem::device_vector<const unsigned int> d_node_eta_index(
        payload.node_eta_index);
    const vecmem::device_vector<const unsigned int> d_node_phi_index(
        payload.node_phi_index);
    vecmem::device_vector<unsigned int> d_phi_cusums(payload.phi_cusums);
    vecmem::device_vector<float4> d_node_params(payload.node_params);
    vecmem::device_vector<float> d_node_phi(payload.node_phi);
    vecmem::device_vector<unsigned int> d_node_index(payload.node_index);
    const vecmem::device_vector<const unsigned int> d_original_sp_idx(
        payload.original_sp_idx);
    const vecmem::device_vector<const float> d_tau_lut(payload.tau_lut);

    const gbts_sort_nodes_params& ap = payload.gbts_sort_nodes_params;

    const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
    const unsigned int blockDimX = thread_id.getBlockDimX();
    const unsigned int gridDimX = thread_id.getGridDimX();

    for (unsigned int globalIndex = globalIdx; globalIndex < payload.nNodes;
         globalIndex += blockDimX * gridDimX) {

        const float4 sp = d_sp_params[globalIndex];

        const float Phi = math::atan2(sp.y, sp.x);
        const float r = math::sqrt(sp.x * sp.x + sp.y * sp.y);
        const float z = sp.z;

        // Default to the full |tau| acceptance for nodes that carry no usable
        // cluster width (sp.w <= 0); the per-edge cuts then rely on these
        // bounds.
        float min_tau = 0.0f;
        float max_tau = ap.maxTau;

        if (sp.w > 0) {  // type 0 only
            if (ap.useTauLUT) {
                // LUT is laid out as [w_bin_edge, min_tau_0, max_tau_0,
                // min_tau_1, max_tau_1] per bin.
                const int tau_bin =
                    5 * static_cast<int>(
                            math::floor(ap.tau_lut_inv_bin * sp.w) - 1.0f);
                if (tau_bin > -1 && tau_bin < static_cast<int>(ap.tauLutSize)) {
                    min_tau =
                        d_tau_lut[static_cast<unsigned int>(tau_bin) + 1u];
                    max_tau =
                        d_tau_lut[static_cast<unsigned int>(tau_bin) + 2u];
                }
                if (max_tau < 0.0f) {
                    max_tau = ap.maxTau;
                }
                if (min_tau < 0.0f) {
                    min_tau = 0.0f;
                }
            } else {
                // linear fit + correction for short clusters
                min_tau = ap.tMin_slope * (sp.w - ap.offset);
                max_tau = ap.tMax_min +
                          ap.tMax_correction / (sp.w + ap.offset) +
                          ap.tMax_slope * (sp.w - ap.offset);
            }
        }

        const unsigned int eta_index = d_node_eta_index[globalIndex];
        const unsigned int histo_bin =
            d_node_phi_index[globalIndex] + payload.nPhiBins * eta_index;

        const unsigned int pos =
            vecmem::device_atomic_ref<unsigned int>(d_phi_cusums[histo_bin])
                .fetch_add(1);

        d_node_params[pos] = float4{min_tau, max_tau, r, z};
        d_node_phi[pos] = Phi;
        d_node_index[pos] = d_original_sp_idx[globalIndex];
    }
}

}  // namespace traccc::device
