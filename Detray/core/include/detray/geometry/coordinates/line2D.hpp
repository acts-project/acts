// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"

namespace detray {

/// Projection into a line coordinate frame
template <concepts::algebra algebra_t>
struct line2D {
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;

  /// Local point type in line coordinates
  using loc_point = point2_type;

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 3D line point
  DETRAY_HOST_DEVICE
  static inline point3_type global_to_local_3D(const transform3_type &trf,
                                               const point3_type &p,
                                               const vector3_type &dir) {
    const auto local3 = trf.point_to_local(p);

    // Line direction
    const vector3_type z = trf.z();

    // Line center
    const point3_type t = trf.translation();

    // Radial vector
    const vector3_type r = vector::cross(z, dir);

    // Assign the sign depending on the position w.r.t line
    // Right: -1
    // Left: 1
    const scalar_type sign = vector::dot(r, t - p) > 0.f ? -1.f : 1.f;

    return {sign * vector::perp(local3), local3[2], vector::phi(local3)};
  }

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 3D line point
  DETRAY_HOST_DEVICE
  static inline loc_point global_to_local(const transform3_type &trf,
                                          const point3_type &p,
                                          const vector3_type &dir) {
    const point3_type local3 =
        line2D<algebra_t>::global_to_local_3D(trf, p, dir);

    return {local3[0], local3[1]};
  }

  /// This method transforms from a local 3D line point to a point in
  /// the global cartesian 3D frame
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const point3_type &p) {
    const scalar_type R = math::fabs(p[0]);
    const point3_type local = {R * math::cos(p[2]), R * math::sin(p[2]), p[1]};

    return trf.point_to_global(local);
  }

  /// This method transforms from a local 2D line point to a point in
  /// the global cartesian 3D frame
  template <typename mask_t>
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const mask_t & /*mask*/, const loc_point &p,
      const vector3_type &dir) {
    // Line direction
    const vector3_type z = trf.z();

    // Radial vector
    const vector3_type r = vector::cross(z, dir);

    // Local Z position in global cartesian coordinate
    const point3_type locZ_in_global =
        trf.point_to_global(point3_type{scalar_type(0), scalar_type(0), p[1]});

    return locZ_in_global + p[0] * vector::normalize(r);
  }

  /// @returns the normal vector in global coordinates
  template <typename mask_t>
  DETRAY_HOST_DEVICE static inline vector3_type normal(
      const transform3_type &trf, const point2_type & = {},
      const mask_t & = {}) {
    return trf.z();
  }

  /// @returns the normal vector in global coordinates
  DETRAY_HOST_DEVICE static inline vector3_type normal(
      const transform3_type &trf, const point3_type & = {}) {
    return trf.z();
  }
};

}  // namespace detray
