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
#include "traccc/utils/trigonometric_helpers.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

// For each edge, find up to nMaxNei compatible neighbour edges sharing its
// outer node, recording them and marking both edges as "kept". The edge
// parameters are read from the packed short4 buffer ([eta, curv, phi_z,
// phi_w]) and the cuts are evaluated in float.
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_match_graph_edges(
    const thread_id_t& thread_id,
    const gbts_match_graph_edges_payload& payload) {
  const vecmem::device_vector<const short4> d_edge_params(payload.edge_params);
  const vecmem::device_vector<const uint2> d_edge_nodes(payload.edge_nodes);
  const vecmem::device_vector<const unsigned int> d_num_outgoing_edges(
      payload.num_outgoing_edges);
  const vecmem::device_vector<const unsigned int> d_edge_links(
      payload.edge_links);
  vecmem::device_vector<unsigned char> d_num_neighbours(payload.num_neighbours);
  vecmem::device_vector<unsigned int> d_neighbours(payload.neighbours);
  vecmem::device_vector<int> d_reIndexer(payload.reIndexer);

  const float cut_dphi_max = payload.gbts_match_graph_edges_params.cut_dphi_max;
  const float cut_dcurv_max =
      payload.gbts_match_graph_edges_params.cut_dcurv_max;
  const float cut_deta_max = payload.gbts_match_graph_edges_params.cut_deta_max;
  const float cut_ratio_sum_max =
      payload.gbts_match_graph_edges_params.cut_ratio_sum_max;
  const float deta_inflation =
      payload.gbts_match_graph_edges_params.deta_inflation;
  const float less_scattering_curv =
      payload.gbts_match_graph_edges_params.less_scattering_curv;
  const float much_less_scattering_curv =
      payload.gbts_match_graph_edges_params.much_less_scattering_curv;
  const float high_pT_correction =
      payload.gbts_match_graph_edges_params.high_pT_correction;

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  for (unsigned int globalIndex = globalIdx; globalIndex < payload.nEdges;
       globalIndex += blockDimX * gridDimX) {
    const unsigned int sharedNode = d_edge_nodes[globalIndex].x;

    const unsigned int link_begin = d_num_outgoing_edges[sharedNode];
    // the number of edges leaving the sharedNode
    const unsigned int nLinks =
        d_num_outgoing_edges[sharedNode + 1u] - link_begin;
    if (nLinks == 0u) {
      continue;
    }

    const std::pair<float4, bool> params1 =
        payload.edge_params_decoder.decode_edge_params(
            d_edge_params[globalIndex]);
    // [exp_eta, curv, phi_z, phi_w], inflate cuts

    const float eta2 = params1.first.x;
    const float Phi2 = params1.first.z;
    const float curv2 = params1.first.y;

    const unsigned int nei_pos = payload.nMaxNei * globalIndex;

    unsigned char num_nei = 0;

    for (unsigned int k = 0u; k < nLinks;
         k++) {  // loop over potential neighbours

      if (num_nei >= payload.nMaxNei) {
        break;
      }
      const unsigned int edge2_idx = d_edge_links[link_begin + k];

      const std::pair<float4, bool> params2 =
          payload.edge_params_decoder.decode_edge_params(
              d_edge_params[edge2_idx]);

      // adaptive eta cut based on edge length and curvature
      float deta_max =
          cut_deta_max + deta_inflation * (params2.second || params1.second);
      const float curv = 0.5f * math::fabs(curv2 + params2.first.y);
      const float corr = static_cast<float>((curv < less_scattering_curv) +
                                            (curv < much_less_scattering_curv));
      deta_max *= 1.0f - high_pT_correction * corr;
      const float deta_cut_ratio =
          math::fabs(eta2 - params2.first.x) / deta_max;
      if (deta_cut_ratio > 1.0f) {  // bad match
        continue;
      }

      const float dPhi_cut_ratio =
          math::fabs(traccc::detail::wrap_phi(Phi2 - params2.first.w)) /
          cut_dphi_max;
      if (dPhi_cut_ratio > 1.0f) {
        continue;
      }

      const float dcurv_cut_ratio =
          math::fabs(curv2 - params2.first.y) / cut_dcurv_max;
      if (dcurv_cut_ratio > 1.0f) {
        continue;
      }

      if (deta_cut_ratio + dPhi_cut_ratio + dcurv_cut_ratio >
          cut_ratio_sum_max) {
        continue;
      }

      d_neighbours[nei_pos + num_nei] = edge2_idx;
      d_reIndexer[edge2_idx] = 1;
      ++num_nei;
    }

    d_num_neighbours[globalIndex] = num_nei;

    if (num_nei != 0) {
      d_reIndexer[globalIndex] = 1;
      vecmem::device_atomic_ref<unsigned int>(*payload.nConnectionsCounter)
          .fetch_add(static_cast<unsigned int>(num_nei));
    }
  }
}

}  // namespace traccc::device
