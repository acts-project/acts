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
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/cylindrical2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl;

/// A functor to find intersections between straight line and planar surface
template <detray::concepts::soa algebra_t, bool resolve_pos>
struct ray_intersector_impl<cylindrical2D<algebra_t>, algebra_t, resolve_pos> {
  /// Linear algebra types
  /// @{
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;
  /// @}

  template <typename surface_descr_t>
  using intersection_type =
      intersection2D<surface_descr_t, algebra_t, resolve_pos>;

  template <typename other_algebra_t>
  using trajectory_type = detail::ray<other_algebra_t>;

  // Maximum number of solutions this intersector can produce
  static constexpr std::uint8_t n_solutions{2u};

  /// Always includes the intersection position, in order to resolve the mask
  using result_type = darray<
      intersection_point<algebra_t, point3_type, intersection::contains_pos>,
      n_solutions>;

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
  template <typename mask_t, concepts::algebra other_algebra_t>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<other_algebra_t> &ray, const transform3_type &trf,
      const mask_t &mask, const scalar_type /*overstep_tol*/ = 0.f) const {
    // One or both of these solutions might be invalid
    const auto qe = solve_intersection(ray, mask, trf);

    const auto &pos = ray.pos();
    const auto &dir = ray.dir();
    const point3_type ro{pos[0], pos[1], pos[2]};
    const vector3_type rd{dir[0], dir[1], dir[2]};

    result_type results;

    results[0].path = qe.smaller();
    results[0].point = ro + qe.smaller() * rd;

    results[1].path = qe.larger();
    results[1].point = ro + qe.larger() * rd;

    return results;
  }

 protected:
  /// Calculates the distance to the (two) intersection points on the
  /// cylinder in global coordinates.
  ///
  /// @returns a quadratic equation object that contains the solution(s).
  template <typename mask_t, typename other_algebra_t>
  DETRAY_HOST_DEVICE inline auto solve_intersection(
      const trajectory_type<other_algebra_t> &ray, const mask_t &mask,
      const transform3_type &trf) const {
    const vector3_type &sz = trf.z();
    const point3_type &sc = trf.translation();

    const scalar_type r = mask[mask_t::shape::e_r];

    const auto &pos = ray.pos();
    const auto &dir = ray.dir();
    const point3_type ro{pos[0], pos[1], pos[2]};
    const vector3_type rd{dir[0], dir[1], dir[2]};

    const vector3_type tmp = ro - sc;
    const auto pc_cross_sz = vector::cross(tmp, sz);
    const auto rd_cross_sz = vector::cross(rd, sz);
    const scalar_type a = vector::dot(rd_cross_sz, rd_cross_sz);
    const scalar_type b = 2.f * vector::dot(rd_cross_sz, pc_cross_sz);
    const scalar_type c = vector::dot(pc_cross_sz, pc_cross_sz) - (r * r);

    return detail::quadratic_equation<scalar_type>{a, b, c};
  }
};

}  // namespace detray
