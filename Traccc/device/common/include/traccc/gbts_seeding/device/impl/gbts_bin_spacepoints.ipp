/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

// System include(s).
#include <climits>
#include <utility>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_bin_spacepoints(
    const thread_id_t& thread_id, const gbts_bin_spacepoints_payload& payload) {
  const vecmem::device_vector<const float4> reducedSP(payload.reducedSP);
  const vecmem::device_vector<const unsigned short> spacepointsLayer(
      payload.spacepointsLayer);
  vecmem::device_vector<unsigned int> layerCounts(payload.layerCounts);
  vecmem::device_vector<float4> sp_params(payload.sp_params);
  vecmem::device_vector<unsigned int> original_sp_idx(payload.original_sp_idx);
  const vecmem::device_vector<const std::pair<unsigned int, unsigned int>>
      d_layer_info(payload.layer_info);
  const vecmem::device_vector<const std::pair<float, float>> d_layer_geo(
      payload.layer_geo);
  vecmem::device_vector<unsigned int> d_node_eta_index(payload.node_eta_index);
  vecmem::device_vector<unsigned int> d_node_phi_index(payload.node_phi_index);
  vecmem::device_vector<unsigned int> d_eta_phi_histo(payload.eta_phi_histo);

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  for (unsigned int globalIndex = globalIdx; globalIndex < payload.nSp;
       globalIndex += blockDimX * gridDimX) {
    // --- Stage 1: bin_sp_by_layer
    // ------------------------------------------
    const float4 sp = reducedSP[globalIndex];
    if (sp.w < -CHAR_MAX) {
      continue;
    }
    const unsigned short layerIdx = spacepointsLayer[globalIndex];
    const unsigned int binedIdx =
        vecmem::device_atomic_ref<unsigned int>(layerCounts[layerIdx])
            .fetch_sub(1) -
        1u;
    original_sp_idx[binedIdx] = globalIndex;
    sp_params[binedIdx] = sp;

    // --- Stage 2: node_eta_binning -----------
    const unsigned int layerIdx_u = static_cast<unsigned int>(layerIdx);
    const std::pair<unsigned int, unsigned int> layerInfo =
        d_layer_info[layerIdx_u];
    const unsigned int bin0 = layerInfo.first;
    const unsigned int num_eta_bins = layerInfo.second;
    unsigned int eta_index;
    if (num_eta_bins == 1u) {
      eta_index = bin0;
    } else {
      const std::pair<float, float> layerGeo = d_layer_geo[layerIdx_u];
      const float min_eta = layerGeo.first;
      const float eta_bin_width = layerGeo.second;
      const float r = math::sqrt(sp.x * sp.x + sp.y * sp.y);
      const float t1 = sp.z / r;
      const float eta = -math::log(math::sqrt(1.0f + t1 * t1) - t1);
      const unsigned int binIdx = static_cast<unsigned int>(
          math::max(0.0f, math::min((eta - min_eta) / eta_bin_width,
                                    static_cast<float>(num_eta_bins - 1u))));
      eta_index = bin0 + binIdx;
    }
    d_node_eta_index[binedIdx] = eta_index;

    // --- Stage 3: eta_phi_histo --------------
    const float inv_phiSliceWidth =
        1.0f /
        (traccc::device::TWO_PI_F / static_cast<float>(payload.nPhiBins));
    const float Phi = math::atan2(sp.y, sp.x);
    unsigned int phiIdx = static_cast<unsigned int>(
        (Phi + traccc::device::PI_F) * inv_phiSliceWidth);
    if (phiIdx >= payload.nPhiBins) {
      phiIdx -= payload.nPhiBins;
    }
    d_node_phi_index[binedIdx] = phiIdx;

    const unsigned int histo_bin = phiIdx + payload.nPhiBins * eta_index;
    vecmem::device_atomic_ref<unsigned int>(d_eta_phi_histo[histo_bin])
        .fetch_add(1u);
  }
}

}  // namespace traccc::device
