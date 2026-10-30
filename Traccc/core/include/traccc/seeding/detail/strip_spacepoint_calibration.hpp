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
  vector3 inner_cross_outer_half_vector{};
  vector3 separation_cross_outer_half_vector{};
  vector3 separation_cross_inner_half_vector{};
  point3 outer_center{};
  vector3 outer_half_vector{};
  /// Zero for spacepoints that are not formed from a strip pair
  unsigned int valid{0u};
};

/// Declare all strip spacepoint calibration collection types
using strip_spacepoint_calibration_collection_types =
    collection_types<strip_spacepoint_calibration>;

/// Derive the calibration quantities of a strip spacepoint
///
/// @param outer_center the centre of the outer strip
/// @param outer_half_vector half of the outer strip, from its centre
/// @param inner_half_vector half of the inner strip, from its centre
/// @param separation the vector from the inner to the outer strip centre
/// @return the calibration record, invalid if the strips are degenerate
///
TRACCC_HOST_DEVICE inline strip_spacepoint_calibration
make_strip_spacepoint_calibration(const point3& outer_center,
                                  const vector3& outer_half_vector,
                                  const vector3& inner_half_vector,
                                  const vector3& separation) {
  strip_spacepoint_calibration result{};
  result.inner_cross_outer_half_vector =
      vector::cross(inner_half_vector, outer_half_vector);
  result.separation_cross_outer_half_vector =
      vector::cross(separation, outer_half_vector);
  result.separation_cross_inner_half_vector =
      vector::cross(separation, inner_half_vector);
  result.outer_center = outer_center;
  result.outer_half_vector = outer_half_vector;
  // The scale of the crossing calculation vanishes for parallel strips
  result.valid = vector::dot(result.inner_cross_outer_half_vector,
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
  return make_strip_spacepoint_calibration(outer.center, 0.5f * outer.direction,
                                           0.5f * inner.direction,
                                           outer.center - inner.center);
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
    const vector3& direction, const strip_spacepoint_calibration& sp,
    point3& calibrated, const scalar tolerance) {
  const scalar scale = vector::dot(direction, sp.inner_cross_outer_half_vector);

  // position along the inner strip, relative to its half length
  const scalar s_inner =
      vector::dot(direction, sp.separation_cross_outer_half_vector);
  if (math::fabs(s_inner) > math::fabs(scale) * tolerance) {
    return false;
  }

  // position along the outer strip, relative to its half length
  const scalar s_outer =
      vector::dot(direction, sp.separation_cross_inner_half_vector);
  if (math::fabs(s_outer) > math::fabs(scale) * tolerance) {
    return false;
  }

  calibrated = sp.outer_center + (s_outer / scale) * sp.outer_half_vector;
  return true;
}

}  // namespace traccc
