/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/container.hpp"
#include "traccc/seeding/detail/strip_geometry.hpp"

// System include(s).
#include <array>

namespace traccc {

/// Direction independent quantities of a spacepoint formed from a strip pair
///
/// A strip spacepoint is only known to half a strip length along its strips.
/// These quantities allow it to be moved along them to the crossing of a given
/// track direction, without repeating the derivation for every seed candidate.
/// They are computed once per spacepoint, so the spacepoint EDM itself stays
/// free of strip geometry.
///
struct strip_spacepoint_calibration {
  std::array<scalar, 3> inner_cross_outer_half_vector{};
  std::array<scalar, 3> separation_cross_outer_half_vector{};
  std::array<scalar, 3> separation_cross_inner_half_vector{};
  std::array<scalar, 3> outer_center{};
  std::array<scalar, 3> outer_half_vector{};
  /// Zero for spacepoints that are not formed from a strip pair
  unsigned int valid{0u};
};

/// Declare all strip spacepoint calibration collection types
using strip_spacepoint_calibration_collection_types =
    collection_types<strip_spacepoint_calibration>;

namespace details {

TRACCC_HOST_DEVICE inline std::array<scalar, 3> strip_cross(
    const std::array<scalar, 3>& a, const std::array<scalar, 3>& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

TRACCC_HOST_DEVICE inline scalar strip_dot(const std::array<scalar, 3>& a,
                                           const std::array<scalar, 3>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

}  // namespace details

/// Derive the calibration quantities of a strip spacepoint
///
/// @param outer_center the centre of the outer strip
/// @param outer_half_vector half of the outer strip, from its centre
/// @param inner_half_vector half of the inner strip, from its centre
/// @param separation the vector from the inner to the outer strip centre
/// @return the calibration record, invalid if the strips are degenerate
///
TRACCC_HOST_DEVICE inline strip_spacepoint_calibration
make_strip_spacepoint_calibration(
    const std::array<scalar, 3>& outer_center,
    const std::array<scalar, 3>& outer_half_vector,
    const std::array<scalar, 3>& inner_half_vector,
    const std::array<scalar, 3>& separation) {

  strip_spacepoint_calibration result{};
  result.inner_cross_outer_half_vector =
      details::strip_cross(inner_half_vector, outer_half_vector);
  result.separation_cross_outer_half_vector =
      details::strip_cross(separation, outer_half_vector);
  result.separation_cross_inner_half_vector =
      details::strip_cross(separation, inner_half_vector);
  result.outer_center = outer_center;
  result.outer_half_vector = outer_half_vector;
  // The scale of the crossing calculation vanishes for parallel strips
  result.valid = details::strip_dot(result.inner_cross_outer_half_vector,
                                    result.inner_cross_outer_half_vector) > 0.f
                     ? 1u
                     : 0u;
  return result;
}

/// Derive the calibration quantities from the strips of the two measurements
/// forming the spacepoint
///
/// @param outer the strip of the spacepoint's first measurement
/// @param inner the strip of the spacepoint's second measurement
/// @return the calibration record, invalid if either strip is invalid
///
TRACCC_HOST_DEVICE inline strip_spacepoint_calibration
make_strip_spacepoint_calibration(const details::strip_material& outer,
                                  const details::strip_material& inner) {

  if ((outer.valid == 0u) || (inner.valid == 0u)) {
    return {};
  }
  // strip_material::direction spans the full strip
  const std::array<scalar, 3> outer_half_vector{0.5f * outer.direction[0],
                                                0.5f * outer.direction[1],
                                                0.5f * outer.direction[2]};
  const std::array<scalar, 3> inner_half_vector{0.5f * inner.direction[0],
                                                0.5f * inner.direction[1],
                                                0.5f * inner.direction[2]};
  const std::array<scalar, 3> outer_center{outer.center[0], outer.center[1],
                                           outer.center[2]};
  const std::array<scalar, 3> separation{outer.center[0] - inner.center[0],
                                         outer.center[1] - inner.center[1],
                                         outer.center[2] - inner.center[2]};
  return make_strip_spacepoint_calibration(outer_center, outer_half_vector,
                                           inner_half_vector, separation);
}

/// Move a strip spacepoint along its strips to the crossing of a track with
/// the given direction
///
/// @param[in] direction the track direction at the spacepoint
/// @param[in] sp the calibration quantities of the spacepoint
/// @param[out] calibrated the position of the crossing
/// @param[in] tolerance allowed excess over the strip half length
/// @return false if the crossing is not on the strips
///
TRACCC_HOST_DEVICE inline bool calibrate_strip_spacepoint(
    const std::array<scalar, 3>& direction,
    const strip_spacepoint_calibration& sp, std::array<scalar, 3>& calibrated,
    const scalar tolerance) {

  const scalar scale =
      details::strip_dot(direction, sp.inner_cross_outer_half_vector);

  // position along the inner strip, relative to its half length
  const scalar s_inner =
      details::strip_dot(direction, sp.separation_cross_outer_half_vector);
  if (math::fabs(s_inner) > math::fabs(scale) * tolerance) {
    return false;
  }

  // position along the outer strip, relative to its half length
  const scalar s_outer =
      details::strip_dot(direction, sp.separation_cross_inner_half_vector);
  if (math::fabs(s_outer) > math::fabs(scale) * tolerance) {
    return false;
  }

  const scalar s_outer_norm = s_outer / scale;
  calibrated = {sp.outer_center[0] + sp.outer_half_vector[0] * s_outer_norm,
                sp.outer_center[1] + sp.outer_half_vector[1] * s_outer_norm,
                sp.outer_center[2] + sp.outer_half_vector[2] * s_outer_norm};
  return true;
}

}  // namespace traccc
