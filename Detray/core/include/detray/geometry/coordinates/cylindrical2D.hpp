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

/// Projection into a 2D cylindrical coordinate frame
template <concepts::algebra algebra_t>
struct cylindrical2D {
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;

  /// Local point type in 2D cylindrical coordinates
  using loc_point = point2_type;

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 3D cylindrical point
  DETRAY_HOST_DEVICE
  static inline point3_type global_to_local_3D(const transform3_type &trf,
                                               const point3_type &p,
                                               const vector3_type & /*dir*/) {
    const point3_type local3{trf.point_to_local(p)};
    const scalar_type r{vector::perp(local3)};

    return {r * vector::phi(local3), local3[2], r};
  }

  /// This method transforms a point from a global cartesian 3D frame to a
  /// local 2D cylindrical point
  DETRAY_HOST_DEVICE
  static inline loc_point global_to_local(const transform3_type &trf,
                                          const point3_type &p,
                                          const vector3_type & /*dir*/) {
    const point3_type local3{trf.point_to_local(p)};

    return {vector::perp(local3) * vector::phi(local3), local3[2]};
  }

  /// This method transform from a local 3D cylindrical point to a point in
  /// the global cartesian 3D frame
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const point3_type &p) {
    const scalar_type r{p[2]};
    const scalar_type phi{p[0] / r};
    const scalar_type x{r * math::cos(phi)};
    const scalar_type y{r * math::sin(phi)};
    const scalar_type z{p[1]};

    return trf.point_to_global(point3_type{x, y, z});
  }

  /// This method transform from a local 2D cylindrical point to a point in
  /// the global cartesian 3D frame
  template <typename mask_t>
  DETRAY_HOST_DEVICE static inline point3_type local_to_global(
      const transform3_type &trf, const mask_t &mask, const loc_point &p,
      const vector3_type & /*dir*/) {
    return cylindrical2D<algebra_t>::local_to_global(
        trf, {p[0], p[1], mask[mask_t::shape::e_r]});
  }

  /// @returns the normal vector in global coordinates
  template <typename mask_t>
  DETRAY_HOST_DEVICE static inline vector3_type normal(
      const transform3_type &trf, const point2_type &p, const mask_t &mask) {
    const scalar_type phi{p[0] / mask[mask_t::shape::e_r]};
    const vector3_type local_normal{math::cos(phi), math::sin(phi),
                                    scalar_type(0)};

    // normal vector in global coordinate
    return trf.vector_to_global(local_normal);
  }

  /// @returns the normal vector in global coordinates given a local position
  /// @param p
  DETRAY_HOST_DEVICE static inline vector3_type normal(
      const transform3_type &trf, const point3_type &p) {
    const scalar_type phi{p[0] / p[2]};
    const vector3_type local_normal{math::cos(phi), math::sin(phi),
                                    scalar_type(0)};

    // normal vector in global coordinate
    return trf.vector_to_global(local_normal);
  }

};  // struct cylindrical2D

}  // namespace detray
