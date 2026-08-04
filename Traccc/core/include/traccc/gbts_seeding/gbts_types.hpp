/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"  // traccc::constant
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/qualifiers.hpp"  // TRACCC_HOST_DEVICE, TRACCC_ALIGN

// System include(s).
#include <climits>
#include <cmath>

namespace traccc {

// ---------------------------------------------------------------------------
// Minimal fixed-size vector storage types, and constructors for them.
// ---------------------------------------------------------------------------
struct TRACCC_ALIGN(8) float2 {
  float x, y;
};
struct TRACCC_ALIGN(16) float4 {
  float x, y, z, w;
};
struct TRACCC_ALIGN(8) int2 {
  int x, y;
};
struct TRACCC_ALIGN(8) uint2 {
  unsigned int x, y;
};
struct TRACCC_ALIGN(4) short2 {
  short x, y;
};
struct TRACCC_ALIGN(8) short4 {
  short x, y, z, w;
};

namespace device {

inline constexpr float PI_F = traccc::constant<float>::pi;
inline constexpr float TWO_PI_F = 2.0f * traccc::constant<float>::pi;

// for float->short->float edge params conversion
inline constexpr float short_max = static_cast<float>(SHRT_MAX);

struct edge_params_converter {
  float max_curv = 5e-4f;
  float max_eta = 36.0f;

  edge_params_converter(const float maxcurv, const float maxtau) {
    max_curv = maxcurv;
    max_eta = -1.0f * math::log(math::sqrt(1.0f + maxtau * maxtau) - maxtau) *
              1.0001f;
  }

  inline TRACCC_HOST_DEVICE short4 make_edge_params(float eta, float curv,
                                                    float Phi2, float Phi1,
                                                    bool long_edge) const {
    short4 params = short4{static_cast<short>(short_max * eta / max_eta),
                           static_cast<short>(short_max * curv / max_curv),
                           static_cast<short>(short_max * Phi2 / (TWO_PI_F)),
                           static_cast<short>(short_max * Phi1 / (TWO_PI_F))};
    // use last bit of phi1 to signal long_edge
    params.w -= (params.w & 0x1);
    params.w += long_edge;
    return params;
  }

  inline TRACCC_HOST_DEVICE std::pair<float4, bool> decode_edge_params(
      short4 sp) const {
    float4 float_params =
        float4{max_eta * static_cast<float>(sp.x) / short_max,
               max_curv * static_cast<float>(sp.y) / short_max,
               TWO_PI_F * static_cast<float>(sp.z) / short_max,
               TWO_PI_F * static_cast<float>(sp.w) / short_max};
    return std::make_pair(float_params, sp.w & 0x1);
  }
};

}  // namespace device

}  // namespace traccc
