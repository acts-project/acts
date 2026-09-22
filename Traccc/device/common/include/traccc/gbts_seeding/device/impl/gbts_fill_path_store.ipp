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
#include "traccc/gbts_seeding/device/details/gbts_fit_segments.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

// Detray include(s).
#include <detray/utils/find_bound.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fill_path_store(
    const thread_id_t& thread_id, const gbts_fill_path_store_payload& payload) {
  vecmem::device_vector<int2> d_path_store(payload.path_store);
  const vecmem::device_vector<const unsigned int> d_output_graph(
      payload.output_graph);
  const vecmem::device_vector<const unsigned char> d_levels(payload.levels);
  const vecmem::device_vector<const int2> d_outgoing_paths(
      payload.outgoing_paths);
  const vecmem::device_vector<const unsigned int> d_path_counts(
      payload.path_counts);
  vecmem::device_vector<int2> d_seed_proposals(payload.seed_proposals);
  vecmem::device_vector<char> d_seed_ambiguity(payload.seed_ambiguity);
  vecmem::device_vector<unsigned long long int> d_edge_bids(payload.edge_bids);
  const vecmem::device_vector<const float4> d_sp_reduced(payload.reducedSP);
  const gbts_fit_segments_params& fit_params = payload.gbts_fit_segments_params;

  // Row-major output graph: each edge owns a contiguous block of
  // nei_start + max_num_neighbours ints ([node1, node2, nNei, nei0..]).
  const unsigned int edge_size =
      gbts_consts::nei_start + payload.max_num_neighbours;

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  const unsigned int path_count =
      vecmem::device_vector<const unsigned int>(payload.path_count)[0];
  const unsigned int nPaths =
      (path_count < payload.nPathsMax) ? path_count : payload.nPathsMax;
  for (unsigned int path_idx = globalIdx; path_idx < nPaths;
       path_idx += blockDimX * gridDimX) {
    const unsigned int root = static_cast<unsigned int>(
        detray::detail::upper_bound(
            d_path_counts.begin(),
            d_path_counts.begin() + payload.nConnectedEdges, path_idx) -
        d_path_counts.begin());
    const unsigned int base =
        d_path_counts[root] -
        (1u + static_cast<unsigned int>(d_outgoing_paths[root].x));

    // Descend the preorder layout to the edge this path ends at.
    unsigned int cur_edge = root;
    unsigned int cur_path_idx = base;
    unsigned int offset = path_idx - base;
    int parent_path_idx = -1;
    // The edges of the path, root first.
    unsigned int chain[gbts_consts::max_seed_candidate_length + 1u];
    unsigned int depth = 1u;
    chain[0] = root;
    while (offset > 0u) {
      --offset;
      const unsigned int edge_pos = edge_size * cur_edge;
      const unsigned int nNei = d_output_graph[edge_pos + gbts_consts::nNei];
      const unsigned char level = d_levels[cur_edge];
      unsigned int acc = 0u;
      bool found = false;
      for (unsigned int k = 0u; k < nNei; ++k) {
        const unsigned int child =
            d_output_graph[edge_pos + gbts_consts::nei_start + k];
        if (level != d_levels[child] + 1u) {
          continue;
        }
        const unsigned int size =
            1u + static_cast<unsigned int>(d_outgoing_paths[child].x);
        if (offset < size) {
          parent_path_idx = static_cast<int>(cur_path_idx);
          cur_path_idx = cur_path_idx + 1u + acc;
          cur_edge = child;
          if (depth < gbts_consts::max_seed_candidate_length + 1u) {
            chain[depth++] = child;
          }
          found = true;
          break;
        }
        offset -= size;
        acc += size;
      }
      if (!found) {
        break;
      }
    }
    d_path_store[path_idx] = int2{static_cast<int>(cur_edge), parent_path_idx};
    d_seed_proposals[path_idx] = int2{-1, -1};

    // Segment fit of the path, leaf to root. A single edge is never a seed.
    if (parent_path_idx < 0) {
      continue;
    }
    unsigned char length = 1;
    bool toggle = false;
    details::edgeState state1;
    details::edgeState state2;
    const unsigned int leaf_pos = edge_size * chain[depth - 1u];
    const unsigned int nodeidx1 = d_output_graph[leaf_pos + gbts_consts::node1];
    const traccc::float4 node1 = d_sp_reduced[nodeidx1];
    const unsigned int nodeidx2 = d_output_graph[leaf_pos + gbts_consts::node2];
    traccc::float4 node2 = d_sp_reduced[nodeidx2];
    state1.initialize(node2, node1);
    for (unsigned int i = depth - 1u; i > 0u; --i) {
      const unsigned int nodeidx =
          d_output_graph[edge_size * chain[i - 1u] + gbts_consts::node2];
      node2 = d_sp_reduced[nodeidx];
      if (toggle) {
        if (!details::gbts_kalman_update(&state1, &state2, node2, fit_params,
                                         payload.max_z0)) {
          state1 = state2;
          break;
        }
      } else if (!details::gbts_kalman_update(&state2, &state1, node2,
                                              fit_params, payload.max_z0)) {
        break;
      }
      toggle = !toggle;
      length++;
    }
    if (length < payload.minLevel) {
      continue;
    }
    // state1 is the final state
    //  can cut more strongly now the fit is done
    if (math::fabs(state1.m_X[2]) * fit_params.final_curv_cut_tighten *
            fit_params.inv_max_curvature >
        1.0f) {
      continue;
    }
    // prefer seeds that reach to the outer edge of the detector for better
    // resolution at high pT
    if (math::fabs(state1.m_Y[1]) > fit_params.zmax / fit_params.rmax) {
      state1.m_J += fit_params.add_hit * math::fabs(node1.z) / fit_params.zmax;
    } else {
      const float r2_max = (node1.x) * (node1.x) + (node1.y) * (node1.y);
      state1.m_J += fit_params.add_hit * math::sqrt(r2_max) / fit_params.rmax;
    }
    const int qual = static_cast<int>(fit_params.qual_scale * state1.m_J);
    d_seed_proposals[path_idx] = int2{qual, static_cast<int>(path_idx)};

    // Bid for the path's last edge. The loser is marked ambiguous.
    const unsigned long long int seed_bid =
        (static_cast<unsigned long long int>(qual) << 32) |
        static_cast<unsigned long long int>(path_idx);
    const unsigned long long int competing_offer =
        vecmem::device_atomic_ref<unsigned long long int>(d_edge_bids[cur_edge])
            .fetch_max(seed_bid);
    if (competing_offer > seed_bid) {
      d_seed_ambiguity[path_idx] = -1;
    } else if (competing_offer != 0ull) {
      d_seed_ambiguity[static_cast<unsigned int>(competing_offer &
                                                 0xFFFFFFFFull)] = -1;
    }
  }
}

}  // namespace traccc::device
