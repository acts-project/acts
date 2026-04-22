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
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl;

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <detray::concepts::soa algebra_t, bool resolve_pos>
struct ray_intersector_impl<concentric_cylindrical2D<algebra_t>, algebra_t,
                            resolve_pos> {
  /// Linear algebra types
  /// @{
  using algebra_type = algebra_t;
  using value_type = dvalue<algebra_t>;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
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
  static constexpr std::uint8_t n_solutions{1u};

  /// Always includes the intersection position, in order to resolve the mask
  using result_type =
      intersection_point<algebra_t, point2_type, intersection::contains_pos>;

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
      const trajectory_type<other_algebra_t> &ray,
      const transform3_type & /*trf*/, const mask_t &mask,
      const scalar_type overstep_tol = 0.f) const {
    const scalar_type r = mask[mask_t::shape::e_r];

    const auto &pos = ray.pos();
    const auto &dir = ray.dir();

    const value_type rd_perp_2{dir[0] * dir[0] + dir[1] * dir[1]};

    // The ray is parallel to all/any cylinder axes (z-axis)...
    if (rd_perp_2 < std::numeric_limits<value_type>::epsilon()) [[unlikely]] {
      return {};
    }

    // ...otherwise, two solutions should exist, if the descriminator is
    // greater than zero
    const scalar_type ro_perp_2 = pos[0] * pos[0] + pos[1] * pos[1];
    const scalar_type rad_diff = r * r - ro_perp_2;

    const scalar_type rd_perp_inv_2 = 1.f / rd_perp_2;
    const scalar_type k = -rd_perp_inv_2 * (pos[0] * dir[0] + pos[1] * dir[1]);
    const scalar_type discr = rd_perp_inv_2 * rad_diff + k * k;

    // No intersection found for any cylinder
    if (detray::detail::all_of(discr < 0.f)) [[unlikely]] {
      return {};
    }

    const scalar_type sqrt_discr = math::sqrt(discr);
    const scalar_type s1 = k + sqrt_discr;
    const scalar_type s2 = k - sqrt_discr;

    // Take the nearest solution in every lane
    auto is_smaller_sol = math::fabs(s1) < math::fabs(s2);

    scalar_type path = 0.f;
    path(is_smaller_sol) = s1;
    path(!is_smaller_sol) = s2;

    // If any of the the near solutions is outside the overstepping
    // tolerance, take the far solution (if it exists)
    const auto outside_overstep_tol = path < overstep_tol;
    if (detray::detail::any_of(outside_overstep_tol) &&
        detray::detail::all_of((discr > 0.f) || !outside_overstep_tol)) {
      is_smaller_sol =
          (math::fabs(s1) >= math::fabs(s2)) && outside_overstep_tol;
      path(is_smaller_sol) = s1;
      path(!is_smaller_sol) = s2;
    }

    point2_type loc;
    loc[0] = 0.f;
    loc[1] = pos[2] + path * dir[2];

    return {path, loc};
  }
};

}  // namespace detray
