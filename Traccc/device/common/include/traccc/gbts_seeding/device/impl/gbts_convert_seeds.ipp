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
#include "traccc/gbts_seeding/device/details/gbts_seed_conversion_helpers.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

// System include(s).
#include <array>

namespace traccc::device {

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
    unsigned long long int dropout_index = 0;
    char best_for_hit = 0;
    int2 prop = d_seed_proposals[prop_idx];
    if (prop.y < 0) {
      continue;
    }
    if (d_seed_ambiguity[prop_idx] == -2) {
      continue;
    }
    details::Tracklet seed;
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
      const traccc::float2 curv_cot_1 = details::gbts_estimate_seed_params(sps);
      sps[1] = d_sp_params[seed.nodes[(seed.size - 1) / 2]];
      const traccc::float2 curv_cot_2 = details::gbts_estimate_seed_params(sps);
      sps[0] = d_sp_params[seed.nodes[seed.size - 2]];
      const traccc::float2 curv_cot_3 = details::gbts_estimate_seed_params(sps);
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

}  // namespace traccc::device
