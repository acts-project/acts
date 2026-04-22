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

/// Projection into a 3D cylindrical coordinate frame
template <concepts::algebra algebra_t>
struct cylindrical3D {
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;

  /// Local point type in 3D cylindrical coordinates
  using loc_point = point3_type;

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 3D cylindrical point
  DETRAY_HOST_DEVICE
  static inline point3_type global_to_local_3D(const transform3_type &trf,
                                               const point3_type &p,
                                               const vector3_type &dir) {
    return cylindrical3D<algebra_t>::global_to_local(trf, p, dir);
  }

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 3D cylindrical point
  DETRAY_HOST_DEVICE
  static inline loc_point global_to_local(const transform3_type &trf,
                                          const point3_type &p,
                                          const vector3_type & /*dir*/) {
    const auto local3 = trf.point_to_local(p);
    return {vector::perp(local3), vector::phi(local3), local3[2]};
  }

  /// This method transforms from a local 3D cylindrical point to a point in
  /// the global cartesian 3D frame
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const point3_type &p) {
    const scalar_type x{p[0] * math::cos(p[1])};
    const scalar_type y{p[0] * math::sin(p[1])};

    return trf.point_to_global(point3_type{x, y, p[2]});
  }

  /// This method transforms from a local 3D cylindrical point to a point in
  /// the global cartesian 3D frame
  template <typename mask_t>
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const mask_t & /*mask*/, const loc_point &p,
      const vector3_type & /*dir*/) {
    return cylindrical3D<algebra_t>::global_to_local(trf, p);
  }

  /// @returns the normal vector in global coordinates given a local position
  /// @param p
  DETRAY_HOST_DEVICE static inline vector3_type normal(
      const transform3_type &trf, const point3_type &p) {
    const vector3_type local_normal{math::cos(p[1]), math::sin(p[1]), 0.f};

    // normal vector in global coordinate
    return trf.vector_to_global(local_normal);
  }

};  // struct cylindrical3D

}  // namespace detray
