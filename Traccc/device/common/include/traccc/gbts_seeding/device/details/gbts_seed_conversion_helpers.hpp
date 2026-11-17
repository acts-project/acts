// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device::details {

struct Tracklet {
  unsigned int
      nodes[traccc::device::gbts_consts::max_seed_candidate_length + 1];
  int size;
};

template <typename T>
TRACCC_HOST_DEVICE inline void set_seed(traccc::edm::seed<T>& seed,
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

}  // namespace traccc::device::details
