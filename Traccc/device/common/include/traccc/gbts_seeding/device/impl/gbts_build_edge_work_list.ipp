/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

template <concepts::thread_id1 thread_id_t, concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void gbts_build_edge_work_list(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const gbts_build_edge_work_list_payload& payload,
    const gbts_build_edge_work_list_shared_payload& shared_payload) {
  const vecmem::device_vector<const uint2> d_bin_pairs(payload.bin_pairs);
  const vecmem::device_vector<const unsigned int> d_eta_bin_offsets(
      payload.eta_bin_offsets);
  vecmem::device_vector<unsigned int> d_pair_work_begin(
      payload.pair_work_begin);
  const vecmem::device_vector<const float> d_bin_rads(payload.bin_rads);
  vecmem::device_vector<gbts_graph_edges_work_item> d_work_items(
      payload.work_items);
  vecmem::device_vector<unsigned int> scratch(shared_payload.scratch);

  const unsigned int threadIndex = thread_id.getLocalThreadIdX();
  const unsigned int blockSize = thread_id.getBlockDimX();

  // Every thread owns one contiguous strip of the bin pairs; a block scan of
  // the strip totals gives the offsets of the strips.
  const unsigned int strip = (payload.nBinPairs + blockSize - 1u) / blockSize;
  const unsigned int begin = threadIndex * strip;
  const unsigned int end = math::min(begin + strip, payload.nBinPairs);
  // Count the work items in this thread's strip
  unsigned int total = 0u;
  for (unsigned int pair = begin; pair < end; pair++) {
    const uint2 bins = d_bin_pairs[pair];
    const unsigned int n1 =
        d_eta_bin_offsets[bins.x + 1u] - d_eta_bin_offsets[bins.x];
    const unsigned int n2 =
        d_eta_bin_offsets[bins.y + 1u] - d_eta_bin_offsets[bins.y];
    if ((n1 > 0u) && (n2 > 0u)) {
      total += 1u + (n1 - 1u) / gbts_consts::edge_chunk_size;
    }
  }

  // The inclusive block scan of the totals gives the starting index of each
  // thread's strip in the work-item array.
  const unsigned int value = total;
  scratch[threadIndex] = value;
  barrier.blockBarrier();
  for (unsigned int offset = 1u; offset < blockSize; offset *= 2u) {
    const unsigned int add =
        (threadIndex >= offset) ? scratch[threadIndex - offset] : 0u;
    barrier.blockBarrier();
    scratch[threadIndex] += add;
    barrier.blockBarrier();
  }

  unsigned int running = scratch[threadIndex] - total;
  // A local copy: std::min takes references, and the namespace-scope
  // constant cannot be referenced from device code.
  constexpr unsigned int chunk_size = gbts_consts::edge_chunk_size;
  for (unsigned int pair = begin; pair < end; pair++) {
    const uint2 bins = d_bin_pairs[pair];
    const unsigned int inner_begin = d_eta_bin_offsets[bins.x];
    const unsigned int inner_end = d_eta_bin_offsets[bins.x + 1u];
    const unsigned int outer_begin = d_eta_bin_offsets[bins.y];
    const unsigned int outer_end = d_eta_bin_offsets[bins.y + 1u];
    d_pair_work_begin[pair] = running;
    // Skip empty bins
    if ((inner_end == inner_begin) || (outer_end == outer_begin)) {
      continue;
    }
    const float max_delta_r =
        math::fabs(d_bin_rads[2u * bins.y + 1u] - d_bin_rads[2u * bins.x]);

    const gbts_dphi_window_params& dp = payload.gbts_dphi_window_params;
    const float delta_phi =
        (max_delta_r < dp.low_dr_threshold)
            ? dp.min_delta_phi_low_dr + dp.dphi_coeff_low_dr * max_delta_r
            : dp.min_delta_phi + dp.dphi_coeff * max_delta_r;

    const unsigned int chunks =
        1u + (inner_end - inner_begin - 1u) / chunk_size;
    for (unsigned int chunk_idx = 0u; chunk_idx < chunks; chunk_idx++) {
      const unsigned int chunk_begin = inner_begin + chunk_idx * chunk_size;
      d_work_items[running + chunk_idx] = gbts_graph_edges_work_item{
          pair,         // Index of the bin pair
          chunk_idx,    // Index of the chunk within the inner bin
          chunk_begin,  // Index of the first inner node in the chunk
          math::min(
              chunk_size,
              inner_end - chunk_begin),  // Number of inner nodes in the chunk
          outer_begin,  // Index of the first outer node in the outer bin
          outer_end,    // Index of the last outer node in the outer bin
          delta_phi     // Delta-phi window for the bin pair
      };
    }
    running += chunks;
  }
  if (threadIndex == 0u) {
    d_pair_work_begin[payload.nBinPairs] = scratch[blockSize - 1u];
  }
}

}  // namespace traccc::device
