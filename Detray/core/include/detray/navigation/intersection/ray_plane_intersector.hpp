// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
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
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl;

/// A functor to find intersections between straight line and planar surface
template <detray::concepts::aos algebra_t, bool resolve_pos>
struct ray_intersector_impl<cartesian2D<algebra_t>, algebra_t, resolve_pos> {
  using algebra_type = algebra_t;
  using frame_type = cartesian2D<algebra_t>;
  using point_type = dpoint3D<algebra_t>;

  template <typename surface_descr_t>
  using intersection_type =
      intersection2D<surface_descr_t, algebra_t, resolve_pos>;

  template <typename other_algebra_t>
  using trajectory_type = detail::ray<other_algebra_t>;

  // Maximum number of solutions this intersector can produce
  static constexpr std::uint8_t n_solutions{1u};

  /// Always includes the intersection position, in order to resolve the mask
  using result_type =
      intersection_point<algebra_t, point_type, intersection::contains_pos>;

  /// Operator function to find intersections between ray and planar mask
  ///
  /// @param ray is the input ray trajectory
  /// @param sf the surface handle the mask is associated with
  /// @param mask is the input mask that defines the surface extent
  /// @param trf is the surface placement transform
  /// @param mask_tolerance is the tolerance for mask edges
  /// @param overstep_tol negative cutoff for the path
  ///
  /// @return the intersection
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<algebra_t> &ray, const dtransform3D<algebra_t> &trf,
      const dscalar<algebra_t> /*overstep_tol*/ = 0.f) const {
    using scalar_t = dscalar<algebra_t>;
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    // Retrieve the surface normal & translation (context resolved)
    const vector3_t &sn = trf.z();
    const vector3_t &st = trf.translation();

    // Intersection code
    const point3_t &ro = ray.pos();
    const vector3_t &rd = ray.dir();
    const scalar_t denom = vector::dot(rd, sn);

    // this is dangerous
    if (denom == 0.f) [[unlikely]] {
      return {};
    }

    scalar_t s{vector::dot(sn, st - ro) / denom};
    point3_t glob_pos = ro + s * rd;

    return {s, glob_pos};
  }
};

template <detray::concepts::aos algebra_t, bool resolve_pos>
struct ray_intersector_impl<polar2D<algebra_t>, algebra_t, resolve_pos>
    : public ray_intersector_impl<cartesian2D<algebra_t>, algebra_t,
                                  resolve_pos> {};

}  // namespace detray
