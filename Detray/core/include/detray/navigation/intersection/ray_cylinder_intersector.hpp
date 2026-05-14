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
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/quadratic_equation.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl;

/// A functor to find intersections between a ray and a 2D cylinder mask
template <detray::concepts::aos algebra_t, bool resolve_pos>
struct ray_intersector_impl<cylindrical2D<algebra_t>, algebra_t, resolve_pos> {
  using algebra_type = algebra_t;
  using frame_type = cylindrical2D<algebra_t>;
  using point_type = dpoint3D<algebra_t>;

  template <typename surface_descr_t>
  using intersection_type =
      intersection2D<surface_descr_t, algebra_t, resolve_pos>;

  template <typename other_algebra_t>
  using trajectory_type = detail::ray<other_algebra_t>;

  // Maximum number of solutions this intersector can produce
  static constexpr std::uint8_t n_solutions{2u};

  using result_type = darray<
      intersection_point<algebra_t, point_type, intersection::contains_pos>,
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
  template <typename mask_t>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<algebra_t> &ray, const dtransform3D<algebra_t> &trf,
      const mask_t &mask,
      const dscalar<algebra_t> /*overstep_tol*/ = 0.f) const {
    // One or both of these solutions might be invalid
    const auto qe = solve_intersection(ray, mask, trf);

    result_type results;

    switch (qe.solutions()) {
      [[likely]] case 2: {
        results[1].path = qe.larger();
        results[1].point = ray.pos() + qe.larger() * ray.dir();
        // If there are two solutions, reuse the case for a single
        // solution to setup the intersection point with the smaller
        // path in results[0]
        [[fallthrough]];
      }
      case 1: {
        results[0].path = qe.smaller();
        results[0].point = ray.pos() + qe.smaller() * ray.dir();
        break;
      }
      [[unlikely]] default: { /* Do nothing */
      }
    }

    return results;
  }

 protected:
  /// Calculates the distance to the (two) intersection points on the
  /// cylinder in global coordinates.
  ///
  /// @returns a quadratic equation object that contains the solution(s).
  template <typename mask_t>
  DETRAY_HOST_DEVICE inline detail::quadratic_equation<dscalar<algebra_t>>
  solve_intersection(const trajectory_type<algebra_t> &ray, const mask_t &mask,
                     const dtransform3D<algebra_t> &trf) const {
    using scalar_t = dscalar<algebra_t>;
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    const scalar_t r{mask[mask_t::shape::e_r]};
    const vector3_t &sz = trf.z();
    const vector3_t &sc = trf.translation();

    const point3_t &ro = ray.pos();
    const vector3_t &rd = ray.dir();

    const auto pc_cross_sz = vector::cross(ro - sc, sz);
    const auto rd_cross_sz = vector::cross(rd, sz);
    const scalar_t a{vector::dot(rd_cross_sz, rd_cross_sz)};
    const scalar_t b{2.f * vector::dot(rd_cross_sz, pc_cross_sz)};
    const scalar_t c{vector::dot(pc_cross_sz, pc_cross_sz) - (r * r)};

    return detail::quadratic_equation<scalar_t>{a, b, c};
  }
};

}  // namespace detray
