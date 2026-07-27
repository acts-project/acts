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
// parameters are read from the packed float4 buffer ([exp_eta, curv, phi_z,
// phi_w]) and the cuts are evaluated in float.
template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_match_graph_edges(
    const thread_id_t& thread_id,
    const gbts_match_graph_edges_payload& payload) {
  const vecmem::device_vector<const float4> d_edge_params(payload.edge_params);
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
  const float cut_tau_ratio_max =
      payload.gbts_match_graph_edges_params.cut_tau_ratio_max;

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

    const float4 params1 = d_edge_params[globalIndex];  // [exp_eta, curv,
                                                        //  phi_z, phi_w]
    const float uat_2 = 1.0f / params1.x;
    const float Phi2 = params1.z;
    const float curv2 = params1.y;

    const unsigned int nei_pos = payload.nMaxNei * globalIndex;

    unsigned char num_nei = 0;

    for (unsigned int k = 0u; k < nLinks;
         k++) {  // loop over potential neighbours

      if (num_nei >= payload.nMaxNei) {
        break;
      }
      const unsigned int edge2_idx = d_edge_links[link_begin + k];

      const float4 params2 = d_edge_params[edge2_idx];

      const float tau_ratio = params2.x * uat_2 - 1.0f;
      if (math::fabs(tau_ratio) > cut_tau_ratio_max) {  // bad match
        continue;
      }

      const float dPhi = traccc::detail::wrap_phi(Phi2 - params2.w);
      if (math::fabs(dPhi) > cut_dphi_max) {
        continue;
      }

      const float dcurv = curv2 - params2.y;
      if (math::fabs(dcurv) > cut_dcurv_max) {
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
