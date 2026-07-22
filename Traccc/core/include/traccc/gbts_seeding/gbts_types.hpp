/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"      // traccc::constant
#include "traccc/definitions/qualifiers.hpp"  // TRACCC_HOST_DEVICE, TRACCC_ALIGN

// System include(s).
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

namespace device {

inline constexpr float PI_F = traccc::constant<float>::pi;
inline constexpr float TWO_PI_F = 2.0f * traccc::constant<float>::pi;

}  // namespace device

}  // namespace traccc
