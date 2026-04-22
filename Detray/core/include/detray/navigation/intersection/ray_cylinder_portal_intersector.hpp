// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

// System include(s)
#include <type_traits>

namespace detray {

template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl;

/// @brief A functor to find intersections between a straight line and a
/// cylindrical portal surface.
///
/// With the way the navigation works, only the closest one of the two possible
/// intersection points is needed in the case of a cylinderical portal surface.
template <detray::concepts::aos algebra_t, bool resolve_pos>
struct ray_intersector_impl<concentric_cylindrical2D<algebra_t>, algebra_t,
                            resolve_pos> {
  using algebra_type = algebra_t;
  using frame_type = concentric_cylindrical2D<algebra_t>;
  using point_type = dpoint2D<algebra_t>;

  template <typename surface_descr_t>
  using intersection_type =
      intersection2D<surface_descr_t, algebra_t, resolve_pos>;

  template <typename other_algebra_t>
  using trajectory_type = detail::ray<other_algebra_t>;

  // Maximum number of solutions this intersector can produce
  static constexpr std::uint8_t n_solutions{1u};

  using result_type =
      intersection_point<algebra_t, point_type, intersection::contains_pos>;

  /// Operator function to find intersections between ray and cylinder mask
  ///
  /// Intersecting the cylinder from the inside yields one intersection
  /// along the direction of the track and one behind it. These intersections
  /// can be calculated in a simplified way, since the cylinder cannot be
  /// shifted or rotated. Solve: perp(ro + t * rd) = r_cyl
  ///
  /// @param ray is the input ray trajectory
  /// @param sf the surface handle the mask is associated with
  /// @param trf is the surface placement transform
  /// @param mask_tolerance is the tolerance for mask edges
  /// @param overstep_tol negative cutoff for the path
  ///
  /// @return the closest intersection
  template <typename mask_t>
  DETRAY_HOST_DEVICE constexpr result_type point_of_intersection(
      const trajectory_type<algebra_t> &ray,
      const dtransform3D<algebra_t> & /*trf*/, const mask_t &mask,
      const dscalar<algebra_t> overstep_tol = 0.f) const {
    using scalar_t = dscalar<algebra_t>;
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    const scalar_t r{mask[mask_t::shape::e_r]};
    constexpr scalar_t inv{detail::invalid_value<dvalue<algebra_t>>()};

    const point3_t &ro = ray.pos();
    const vector3_t &rd = ray.dir();
    scalar_t path{inv};

    const scalar_t rd_perp_2{rd[0] * rd[0] + rd[1] * rd[1]};

    // The ray is parallel to the cylinder axis (z-axis)...
    if (rd_perp_2 < std::numeric_limits<scalar_t>::epsilon()) [[unlikely]] {
      return {};
    }

    // ...otherwise, two solutions should exist, if the descriminator is
    // greater than zero
    const scalar_t ro_perp_2{ro[0] * ro[0] + ro[1] * ro[1]};

    // Only calculate the path, when not already on surface
    // 'rad_diff' leads to numerical instabilites in single precision on GPU
    if (const double rad_diff{r * r - ro_perp_2};
        math::fabs(rad_diff) > 1.f * unit<scalar_t>::um) {
      const scalar_t rd_perp_inv_2{1.f / rd_perp_2};
      const scalar_t k{-rd_perp_inv_2 * (ro[0] * rd[0] + ro[1] * rd[1])};
      const auto discr{static_cast<scalar_t>(rd_perp_inv_2 * rad_diff + k * k)};

      // No intersection found
      if (discr < 0.f) [[unlikely]] {
        return {};
      }

      const scalar_t sqrt_discr{math::sqrt(discr)};
      const scalar_t s1{k + sqrt_discr};
      const scalar_t s2{k - sqrt_discr};

      // Take the nearest solution
      path = (math::fabs(s1) < math::fabs(s2)) ? s1 : s2;

      // If the near solution is outside the overstepping tolerance, take
      // the far solution (if it exists)
      if (path < overstep_tol && sqrt_discr > 0.f) {
        path = (math::fabs(s1) >= math::fabs(s2)) ? s1 : s2;
      }
    } else {
      path = static_cast<scalar_t>(rad_diff);
    }

    if (path < overstep_tol) {
      path = inv;
    }

    // Only need the global z-component for the mask check
    return {path, point_type{inv, ro[2] + path * rd[2]}};
  }
};

}  // namespace detray
