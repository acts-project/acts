// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"

namespace detray {

/// Projection into a polar coordinate frame
template <concepts::algebra algebra_t>
struct polar2D {
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;

  /// Local point type in polar coordinates
  using loc_point = point2_type;

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 3D polar point
  DETRAY_HOST_DEVICE
  static inline point3_type global_to_local_3D(const transform3_type &trf,
                                               const point3_type &p,
                                               const vector3_type & /*dir*/) {
    const auto local3 = trf.point_to_local(p);
    return {vector::perp(local3), vector::phi(local3), local3[2]};
  }

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 2D polar point
  DETRAY_HOST_DEVICE
  static inline loc_point global_to_local(const transform3_type &trf,
                                          const point3_type &p,
                                          const vector3_type & /*d*/) {
    const auto local3 = trf.point_to_local(p);
    return {vector::perp(local3), vector::phi(local3)};
  }

  /// This method transforms from a local 3D polar point to a point in
  /// the global cartesian 3D frame
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const point3_type &p) {
    const scalar_type x = p[0] * math::cos(p[1]);
    const scalar_type y = p[0] * math::sin(p[1]);

    return trf.point_to_global(point3_type{x, y, p[2]});
  }

  /// This method transforms from a local 2D polar point to a point in
  /// the global cartesian 3D frame
  template <typename mask_t>
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const mask_t & /*mask*/, const loc_point &p,
      const vector3_type & /*dir*/) {
    return polar2D<algebra_t>::local_to_global(trf,
                                               {p[0], p[1], scalar_type(0)});
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
