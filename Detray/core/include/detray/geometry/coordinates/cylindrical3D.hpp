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
  /// Algebra backend used by this coordinate frame.
  using algebra_type = algebra_t;
  /// Scalar type associated with the algebra backend.
  using scalar_type = dscalar<algebra_t>;
  /// Two-dimensional point type.
  using point2_type = dpoint2D<algebra_t>;
  /// Three-dimensional point type.
  using point3_type = dpoint3D<algebra_t>;
  /// Two-dimensional vector type.
  using vector2_type = dvector2D<algebra_t>;
  /// Three-dimensional vector type.
  using vector3_type = dvector3D<algebra_t>;
  /// Three-dimensional transform type.
  using transform3_type = dtransform3D<algebra_t>;

  /// Local point type in 3D cylindrical coordinates
  using loc_point = point3_type;

  /// Convert a global Cartesian point to cylindrical coordinates in the
  /// surface frame.
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Point in global coordinates to transform.
  /// @param dir Direction passed to the coordinate conversion; some frames ignore it.
  /// @returns The point expressed in this local coordinate frame.
  DETRAY_HOST_DEVICE
  static inline point3_type global_to_local_3D(const transform3_type &trf,
                                               const point3_type &p,
                                               const vector3_type &dir) {
    return cylindrical3D<algebra_t>::global_to_local(trf, p, dir);
  }

  /// Convert a global Cartesian point to cylindrical coordinates in the
  /// surface frame.
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Point in global coordinates to transform.
  /// @returns The point expressed in this local coordinate frame.
  DETRAY_HOST_DEVICE
  static inline loc_point global_to_local(const transform3_type &trf,
                                          const point3_type &p,
                                          const vector3_type & /*dir*/) {
    const auto local3 = trf.point_to_local(p);
    return {vector::perp(local3), vector::phi(local3), local3[2]};
  }

  /// Convert cylindrical coordinates in the surface frame to a global
  /// Cartesian point.
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Local cylindrical coordinates (radius, azimuth, and z).
  /// @returns The point expressed in global coordinates.
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const point3_type &p) {
    const scalar_type x{p[0] * math::cos(p[1])};
    const scalar_type y{p[0] * math::sin(p[1])};

    return trf.point_to_global(point3_type{x, y, p[2]});
  }

  /// Apply the mask-aware coordinate conversion for a local point.
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Local cylindrical coordinates.
  /// The mask and direction arguments are accepted to match the generic frame
  /// interface, but this frame does not use them.
  /// @returns The converted point.
  template <typename mask_t>
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const mask_t & /*mask*/, const loc_point &p,
      const vector3_type & /*dir*/) {
    return cylindrical3D<algebra_t>::global_to_local(trf, p);
  }

  /// Evaluate the surface normal at a local cylindrical position.
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Local position where the normal is evaluated.
  /// @returns The surface normal expressed in global coordinates.
  DETRAY_HOST_DEVICE static inline vector3_type normal(
      const transform3_type &trf, const point3_type &p) {
    // normal vector in global coordinate
    return trf.vector_to_global(vector2_type{math::cos(p[1]), math::sin(p[1])});
  }

};  // struct cylindrical3D

}  // namespace detray
