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
#include "detray/geometry/coordinates/concentric_cylindrical2D.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

// System include(s)
#include <tuple>
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
    //
    // rad_diff here is defined as:
    //
    // r * r - (ro[0] * ro[0] + ro[1] * ro[1])
    //
    // But this runs an enormous risk of catastrophic cancellation, which we
    // would like to avoid.
    //
    // To compute this as accurately as possible without using any double
    // precision operations, we'll use two core ideas:
    //
    // 1. The Sterbenz lemma: if x and y are floats of the same sign and
    //    y/2 <= x <= 2y, then x - y is exact and has no rounding error.
    // 2. FMA hardware has no rounding error internally, so if we compute
    //    a = x * y we get a normal, rounded result. If we then compute
    //    e = x * y - a using fma(x, y, -a), we get the exact error of a.
    //
    // We can use 2) to define an exact squaring algorithm, where to square x
    // we compute first a = x * x and then the error e = fma(x, x, -a).
    // Clearly, the exact square of x is then just a + e. We'll call this
    // procedure a, e = fma_square(x).
    //
    // To compute rad_diff, we can thus compute:
    //
    // rsq_a, rsq_e = fma_square(r)
    // ro0sq_a, ro0sq_e = fma_square(ro[0])
    // ro1sq_a, ro1sq_e = fma_square(ro[1])
    //
    // Such that:
    //
    // rad_diff = (rsq_a + rsq_e) - ((ro0sq_a + ro0sq_e) + (ro1sq_a + ro1sq_e))
    //
    // rad_diff = (rsq_a + rsq_e) - (ro0sq_a + ro0sq_e) - (ro1sq_a + ro1sq_e)
    //
    // However, this is not exact by the Sterbenz lemma 1). Full exactness is
    // not possible, but we can still try to minimize error. We know that the
    // magnitudes of the errors will be much smaller than the magnitudes of the
    // approximated squares, so we can first rearrange:
    //
    // rad_diff = (rsq_a - ro0sq_a - ro1sq_a) + (rsq_e - ro0sq_e - ro1sq_e)
    //
    // What matters now is the way we associate the subtractions inside the
    // two sets of brackets. Arguing informally, we are more likely to get
    // exact (or at least less erroneous) results if the scale of the two
    // numbers are as close as possible, so we will first find the bigger
    // and the smaller value of ro:
    //
    // robig = max(ro[0], ro[1])
    // rosmall = min(ro[0], ro[1])
    //
    // Such that we can express the most accurate result possible:
    //
    // rad_diff = ((rsq_a - robigsq_a) - rosmallsq_a) +
    //            ((rsq_e - robigsq_e) - rosmallsq_e)
    const auto fma_square = [](const scalar_t x) {
      scalar_t a = x * x;
      scalar_t e = math::fma(x, x, -a);
      return std::tuple<scalar_t, scalar_t>(a, e);
    };

    const auto [rsq_a, rsq_e] = fma_square(r);
    const auto [robigsq_a, robigsq_e] =
        fma_square(math::max(math::abs(ro[0]), math::abs(ro[1])));
    const auto [rosmlsq_a, rosmlsq_e] =
        fma_square(math::min(math::abs(ro[0]), math::abs(ro[1])));

    const scalar_t rad_diff =
        ((rsq_a - robigsq_a) - rosmlsq_a) + ((rsq_e - robigsq_e) - rosmlsq_e);

    // Only calculate the path, when not already on surface
    if (math::fabs(rad_diff) > 1.f * unit<scalar_t>::um) {
      const scalar_t rd_perp_inv_2{1.f / rd_perp_2};
      const scalar_t k{-rd_perp_inv_2 * (ro[0] * rd[0] + ro[1] * rd[1])};
      const scalar_t c{rd_perp_inv_2 * rad_diff};
      const scalar_t discr{c + k * k};

      // No intersection found
      if (discr < 0.f) [[unlikely]] {
        return {};
      }

      // The options are k +- sqrt(discr). Clearly the square root of the
      // discriminant is always positive. This means that if k is positive,
      // then k + sqrt(discr) will have a greater absolute value than
      // k - sqrt(discr). If k is negative, the logic is reversed.
      const scalar_t sqrt_discr{math::sqrt(discr)};
      const scalar_t sfar{k + math::copysign(sqrt_discr, k)};
      const scalar_t snear{sfar != 0.f ? -c / sfar : sfar};

      path = snear >= overstep_tol ? snear : sfar;
    } else {
      path = rad_diff;
    }

    if (path < overstep_tol) {
      path = inv;
    }

    // Only need the global z-component for the mask check
    return {path, point_type{inv, ro[2] + path * rd[2]}};
  }
};

}  // namespace detray
