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
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl;

/// A functor to find intersections between straight line and planar surface
template <detray::concepts::soa algebra_t, bool resolve_pos>
struct ray_intersector_impl<cartesian2D<algebra_t>, algebra_t, resolve_pos> {
  /// Linear algebra types
  /// @{
  using algebra_type = algebra_t;
  /// Scalar calculation type for the configured algebra layout.
  using scalar_type = dscalar<algebra_t>;
  /// Three-dimensional point type used by this intersector.
  using point3_type = dpoint3D<algebra_t>;
  /// Three-dimensional vector type used by this intersector.
  using vector3_type = dvector3D<algebra_t>;
  /// Transform type used to place the surface.
  using transform3_type = dtransform3D<algebra_t>;
  /// @}

  /// Intersection result type associated with a surface descriptor.
  template <typename surface_descr_t>
  using intersection_type =
      intersection2D<surface_descr_t, algebra_t, resolve_pos>;

  /// Trajectory representation accepted by this intersector.
  template <typename other_algebra_t>
  using trajectory_type = detail::ray<other_algebra_t>;

  // Maximum number of solutions this intersector can produce
  static constexpr std::uint8_t n_solutions{1u};

  /// Always includes the intersection position, in order to resolve the mask
  using result_type =
      intersection_point<algebra_t, point3_type, intersection::contains_pos>;

  /// Operator function to find intersections between ray and planar mask
  ///
  /// @param ray is the input ray trajectory
  /// @param trf is the surface placement transform
  /// @return the intersection
  template <concepts::algebra other_algebra_t>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<other_algebra_t> &ray, const transform3_type &trf,
      const scalar_type /*overstep_tol*/ = 0.f) const {
    // Retrieve the surface normal & translation (context resolved)
    const vector3_type &sn = trf.z();
    const vector3_type &st = trf.translation();

    // Broadcast ray data
    const auto &pos = ray.pos();
    const auto &dir = ray.dir();
    const vector3_type ro{pos[0], pos[1], pos[2]};
    const vector3_type rd{dir[0], dir[1], dir[2]};

    const scalar_type denom = vector::dot(rd, sn);
    const vector3_type diff = st - ro;
    const scalar_type s = vector::dot(sn, diff) / denom;

    const point3_type glob_pos = ro + s * rd;

    return {s, glob_pos};
  }
};

template <detray::concepts::soa algebra_t, bool resolve_pos>
struct ray_intersector_impl<polar2D<algebra_t>, algebra_t, resolve_pos>
    : public ray_intersector_impl<cartesian2D<algebra_t>, algebra_t,
                                  resolve_pos> {};

}  // namespace detray
