// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <array>

namespace traccc::device {

namespace detail {

struct Tracklet {
  unsigned int
      nodes[traccc::device::gbts_consts::max_seed_candidate_length + 1];
  int size;
};

template <typename T>
TRACCC_HOST_DEVICE inline void set_seed(traccc::edm::seed<T> seed,
                                        const int bottom, const int middle,
                                        const int top, const float quality) {
  seed.bottom_index() = bottom;
  seed.middle_index() = middle;
  seed.top_index() = top;
  seed.quality() = quality;
}

TRACCC_HOST_DEVICE inline traccc::float2 gbts_estimate_seed_params(
    const std::array<traccc::float4, 3>& sps) {
  float u[2], v[2];

  const float x0 = sps[1].x;
  const float y0 = sps[1].y;
  const float r0 = math::sqrt(x0 * x0 + y0 * y0);
  const float cosA = x0 / r0;
  const float sinA = y0 / r0;

  for (unsigned int k = 0; k < 2; k++) {
    const unsigned int sp_idx = (k == 1) ? 2u : k;
    const float dx = sps[sp_idx].x - x0;
    const float dy = sps[sp_idx].y - y0;
    const float r2_inv = 1.0f / (dx * dx + dy * dy);
    const float xn = dx * cosA + dy * sinA;
    const float yn = -dx * sinA + dy * cosA;
    u[k] = xn * r2_inv;
    v[k] = yn * r2_inv;
  }

  const float du = u[0] - u[1];
  if (du == 0.0f) {
    return float2{0.0f, 0.0f};
  }
  const float A = (v[0] - v[1]) / du;
  const float B = v[1] - A * u[1];
  const float curv =
      1000.0f * B / math::sqrt(1 + A * A);  // Curvature from mm^-1 to m^-1
  const float cot_t =
      (sps[2].z - sps[1].z) /
      (math::sqrt(sps[2].x * sps[2].x + sps[2].y * sps[2].y) - r0);
  return float2{curv, cot_t};
}

}  // namespace detail

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_convert_seeds(
    const thread_id_t& thread_id, const gbts_convert_seeds_payload& payload) {
  edm::seed_collection::device seeds_device(payload.output_seeds);
  const vecmem::device_vector<const int2> d_seed_proposals(
      payload.seed_proposals);
  vecmem::device_vector<long long unsigned int> d_seed_bids(payload.seed_bids);
  const vecmem::device_vector<const char> d_seed_ambiguity(
      payload.seed_ambiguity);
  const vecmem::device_vector<const int2> d_path_store(payload.path_store);
  const vecmem::device_vector<const uint2> d_output_edge_nodes(
      payload.output_edge_nodes);
  const vecmem::device_vector<const float4> d_sp_params(payload.reducedSP);
  vecmem::device_vector<unsigned long long int> d_hit_bids(payload.hit_bids);

  const float dcurv_cut_m = payload.gbts_convert_seeds_params.dropout_dcurv_m;
  const float force_dropout_max_curv_m =
      payload.gbts_convert_seeds_params.force_dropout_max_curv_m;
  const float best_hit_frac = payload.gbts_convert_seeds_params.best_hit_frac;
  const float tight_bid_cot_threshold =
      payload.gbts_convert_seeds_params.tight_bid_cot_threshold;
  const bool use_dropout = payload.gbts_convert_seeds_params.use_dropout;

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  const unsigned int nPaths = math::min(
      payload.nPathsMax,
      vecmem::device_vector<const unsigned int>(payload.path_count)[0]);

  for (unsigned int prop_idx = globalIdx; prop_idx < nPaths;
       prop_idx += blockDimX * gridDimX) {
    char diff_code = 0;
    bool force_dropout = false;
    unsigned long long int dropout_index = -1;
    char best_for_hit = 0;
    int2 prop = d_seed_proposals[prop_idx];
    if (prop.y < 0) {
      continue;
    }
    if (d_seed_ambiguity[prop_idx] == -2) {
      continue;
    }
    detail::Tracklet seed;
    seed.size = 0;
    int2 path = int2{0, prop.y};
    while (path.y >= 0) {
      path = d_path_store[static_cast<unsigned int>(path.y)];
      seed.nodes[seed.size++] =
          d_output_edge_nodes[static_cast<unsigned int>(path.x)].x;
      best_for_hit +=
          (prop_idx == (d_hit_bids[seed.nodes[seed.size - 1]] & 0xFFFFFFFFLL));
    }
    seed.nodes[seed.size++] =
        d_output_edge_nodes[static_cast<unsigned int>(path.x)].y;
    best_for_hit +=
        (prop_idx == (d_hit_bids[seed.nodes[seed.size - 1]] & 0xFFFFFFFFLL));

    if (best_for_hit < best_hit_frac * static_cast<float>(seed.size)) {
      continue;
    }
    if (use_dropout) {
      std::array<traccc::float4, 3> sps = {
          d_sp_params[seed.nodes[seed.size - 1]],
          d_sp_params[seed.nodes[(seed.size - 1) / 2 + 1]],
          d_sp_params[seed.nodes[0]]};
      const traccc::float2 curv_cot_1 = detail::gbts_estimate_seed_params(sps);
      sps[1] = d_sp_params[seed.nodes[(seed.size - 1) / 2]];
      const traccc::float2 curv_cot_2 = detail::gbts_estimate_seed_params(sps);
      sps[0] = d_sp_params[seed.nodes[seed.size - 2]];
      const traccc::float2 curv_cot_3 = detail::gbts_estimate_seed_params(sps);
      if ((best_for_hit < seed.size - 1) &
          (fabsf(curv_cot_1.y + curv_cot_2.y +
                 curv_cot_3.y) <  // Checking against the average
                                  // cot(theta) of the three tracklets
           3.0f * tight_bid_cot_threshold) &
          (seed.size < 5)) {  // Don't apply dropout to seeds of length 5
                              // or more. To avoid dropping good seeds.
        continue;
      }
      std::array<float, 3> diff = {fabsf(curv_cot_1.x - curv_cot_2.x),
                                   fabsf(curv_cot_2.x - curv_cot_3.x),
                                   fabsf(curv_cot_1.x - curv_cot_3.x)};
      diff_code = static_cast<char>(4 * (diff[0] < dcurv_cut_m) +
                                    2 * (diff[1] < dcurv_cut_m) +
                                    (diff[2] < dcurv_cut_m));
      force_dropout = fabsf(curv_cot_1.x + curv_cot_2.x + curv_cot_3.x) <
                      3.0f * force_dropout_max_curv_m;
      force_dropout |= (fabsf(curv_cot_1.y + curv_cot_2.y + curv_cot_3.y) <
                        3.0f * tight_bid_cot_threshold) &
                       (diff_code == 0);
    }
    // zero out the last two bits of quality for the dropout index
    // dropout index starts at 2 to prioitize the default/longer seed
    const unsigned int quality_int = prop.x & 0xFFFFFFFC;
    float quality = static_cast<float>(quality_int);
    unsigned long long int seed_bid =
        (static_cast<unsigned long long int>(quality_int) << 32) |
        static_cast<unsigned long long int>(prop.y);
    // use one seed from a consistent pair/set + the inconsistent one
    // sample spacepoints from tracklet to create seeds
    // include 1st order unless either 2 or 3 are consistent with the other
    // and 1
    if (((diff_code != 3) & (diff_code != 6)) | force_dropout) {
      if (payload.gbts_convert_seeds_params.do_seed_sorting) {
        dropout_index = 2;
        d_seed_bids.push_back(seed_bid + (dropout_index << 32));
      } else {
        seeds_device.push_back({seed.nodes[seed.size - 1],
                                seed.nodes[(seed.size - 1) / 2 + 1],
                                seed.nodes[0], quality});
      }
    }
    // include 2nd order if it consistent with 1 and 3 or only 1 and 3 are
    // consistent
    if ((diff_code == 1) | (diff_code == 6)) {
      if (payload.gbts_convert_seeds_params.do_seed_sorting) {
        dropout_index = 1;
        d_seed_bids.push_back(seed_bid + (dropout_index << 32));
      } else {
        seeds_device.push_back({seed.nodes[seed.size - 1],
                                seed.nodes[(seed.size - 1) / 2], seed.nodes[0],
                                quality});
      }
    }
    // include 3rd order if it is consistent with 1 and 2 or only 1 and 2
    // are consistent or if only 2 and 3 are consistent
    if ((diff_code == 2) | (diff_code == 3) | (diff_code == 4) |
        force_dropout) {
      if (payload.gbts_convert_seeds_params.do_seed_sorting) {
        dropout_index = 0;
        d_seed_bids.push_back(seed_bid + (dropout_index << 32));
      } else {
        seeds_device.push_back({seed.nodes[seed.size - 2],
                                seed.nodes[(seed.size - 1) / 2], seed.nodes[0],
                                quality});
      }
    }
  }
}

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fill_seeds_from_sorted_bids(
    const thread_id_t& thread_id,
    const gbts_fill_seeds_from_sorted_bids_payload& payload) {
  edm::seed_collection::device seeds_device(payload.output_seeds);
  vecmem::device_vector<long long unsigned int> d_seed_bids(payload.seed_bids);
  const vecmem::device_vector<const int2> d_path_store(payload.path_store);
  const vecmem::device_vector<const uint2> d_output_edge_nodes(
      payload.output_edge_nodes);

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  for (unsigned int prop_idx = globalIdx; prop_idx < payload.nSeeds;
       prop_idx += blockDimX * gridDimX) {
    // take the dropout index from the saved bid
    unsigned long long int bid = d_seed_bids.at(prop_idx);
    int2 prop = {bid >> 32, bid & 0xFFFFFFFF};
    unsigned int dropout_index = prop.x & 0x3;
    // match quality to non-sorted case
    prop.x -= dropout_index;

    detail::Tracklet seed;
    seed.size = 0;
    int2 path = int2{0, prop.y};
    while (path.y >= 0) {
      path = d_path_store[static_cast<unsigned int>(path.y)];
      seed.nodes[seed.size++] =
          d_output_edge_nodes[static_cast<unsigned int>(path.x)].x;
    }
    seed.nodes[seed.size++] =
        d_output_edge_nodes[static_cast<unsigned int>(path.x)].y;

    const int node1 = seed.nodes[seed.size - 1 - (dropout_index == 0)];
    const int node2 = seed.nodes[(seed.size - 1) / 2 + (dropout_index == 2)];
    const int node3 = seed.nodes[0];

    auto out = seeds_device.at(prop_idx);
    detail::set_seed(out, node1, node2, node3, static_cast<float>(prop.x));
  }
}

}  // namespace traccc::device
