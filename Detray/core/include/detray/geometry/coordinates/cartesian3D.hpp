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

namespace detray {

/// Projection into a 3D cartesian coordinate frame
template <concepts::algebra algebra_t>
struct cartesian3D {
  /// Algebra backend used by this coordinate frame.
  using algebra_type = algebra_t;
  /// Scalar type associated with the algebra backend.
  using scalar_type = dscalar<algebra_t>;
  /// Two-dimensional point type.
  using point2_type = dpoint2D<algebra_t>;
  /// Three-dimensional point type.
  using point3_type = dpoint3D<algebra_t>;
  /// Three-dimensional vector type.
  using vector3_type = dvector3D<algebra_t>;
  /// Three-dimensional transform type.
  using transform3_type = dtransform3D<algebra_t>;

  /// Local point type in 3D cartesian coordinates
  using loc_point = point3_type;

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 3D cartesian point
  /// Convert a point from the global frame to this local coordinate frame.
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Point in global coordinates to transform.
  /// @param dir Direction passed to the coordinate conversion; some frames ignore it.
  /// @returns The point expressed in this local coordinate frame.
  DETRAY_HOST_DEVICE
  static inline point3_type global_to_local_3D(const transform3_type &trf,
                                               const point3_type &p,
                                               const vector3_type &dir) {
    return cartesian3D<algebra_t>::global_to_local(trf, p, dir);
  }

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 3D cartesian point
  /// Convert a point from the global frame to this local coordinate frame.
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Point in global coordinates to transform.
  /// @returns The point expressed in this local coordinate frame.
  DETRAY_HOST_DEVICE
  static inline loc_point global_to_local(const transform3_type &trf,
                                          const point3_type &p,
                                          const vector3_type & /*dir*/) {
    return trf.point_to_local(p);
  }

  /// This method transforms from a local 3D cartesian point to a point in
  /// the global cartesian 3D frame
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Point in local coordinates to transform.
  /// @returns The point expressed in global coordinates.
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const point3_type &p) {
    return trf.point_to_global(p);
  }

  /// This method transforms from a local 3D cartesian point to a point in
  /// the global cartesian 3D frame
  /// @param trf Surface transform defining the local frame in global coordinates.
  /// @param p Point in local coordinates to transform.
  /// @returns The point expressed in global coordinates.
  template <typename mask_t>
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const mask_t & /*mask*/, const loc_point &p,
      const vector3_type & /*dir*/) {
    return cartesian3D<algebra_t>::local_to_global(trf, p);
  }

  /// @returns the normal vector in global coordinates given a local position
  /// @returns The surface normal expressed in global coordinates.
  DETRAY_HOST_DEVICE static inline vector3_type normal(
      const transform3_type & /*trf*/, const point3_type & /*p*/) {
    // Not implemented!
    assert(false);
    return {};
  }

};  // struct cartesian3D

}  // namespace detray
