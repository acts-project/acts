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
#include "detray/geometry/coordinates/line2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl;

/// A functor to find intersections between straight line and planar surface
template <detray::concepts::soa algebra_t, bool resolve_pos>
struct ray_intersector_impl<line2D<algebra_t>, algebra_t, resolve_pos> {
  /// Linear algebra types
  /// @{
  using algebra_type = algebra_t;
  /// Scalar calculation type for the configured algebra layout.
  using scalar_type = dscalar<algebra_t>;
  /// Three-dimensional point type for the returned closest approach.
  using point3_type = dpoint3D<algebra_t>;
  /// Three-dimensional vector type for ray directions.
  using vector3_type = dvector3D<algebra_t>;
  /// Transform type defining the line placement.
  using transform3_type = dtransform3D<algebra_t>;
  /// @}

  /// Intersection result type associated with a surface descriptor.
  template <typename surface_descr_t>
  using intersection_type =
      intersection2D<surface_descr_t, algebra_t, resolve_pos>;

  /// Ray trajectory type accepted by this intersector.
  template <typename other_algebra_t>
  using trajectory_type = detail::ray<other_algebra_t>;

  /// This closest-approach calculation produces one solution.
  static constexpr std::uint8_t n_solutions{1u};

  /// Always includes the intersection position, in order to resolve the mask
  using result_type =
      intersection_point<algebra_t, point3_type, intersection::contains_pos>;

  /// Find the point of closest approach between a ray and the line axis.
  /// @param ray Ray trajectory to compare with the line.
  /// @param trf Transform whose z axis and translation define the line.
  /// The overstep tolerance argument is accepted for interface compatibility
  /// but is unused by the closest-approach calculation.
  /// @returns The path and global point of closest approach.
  template <concepts::algebra other_algebra_t>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<other_algebra_t> &ray, const transform3_type &trf,
      const scalar_type /*overstep_tol*/ = 0.f) const {
    // line direction
    const vector3_type &sz = trf.z();

    // line center
    const point3_type &st = trf.translation();

    // Broadcast ray data
    const auto &pos = ray.pos();
    const auto &dir = ray.dir();
    const vector3_type ro{pos[0], pos[1], pos[2]};
    const vector3_type rd{dir[0], dir[1], dir[2]};

    // Projection of line to track direction
    const scalar_type zd = vector::dot(sz, rd);

    const scalar_type denom = 1.f - (zd * zd);

    // Case for wire is parallel to track
    if (detray::detail::all_of(denom < 1e-5f)) {
      return {};
    }

    // vector from track position to line center
    const vector3_type t2l = st - ro;

    // t2l projection on line direction
    const scalar_type t2l_on_line = vector::dot(t2l, sz);

    // t2l projection on track direction
    const scalar_type t2l_on_track = vector::dot(t2l, rd);

    // path length to the point of closest approach on the track
    scalar_type s = (t2l_on_track - t2l_on_line * zd) / denom;

    // point of closest approach on the track
    point3_type glob_pos = ro + s * rd;

    return {s, glob_pos};
  }
};

}  // namespace detray
