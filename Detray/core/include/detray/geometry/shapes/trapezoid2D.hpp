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
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/cartesian3D.hpp"
#include "detray/geometry/detail/shape_utils.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string_view>

namespace detray {

/// @brief Geometrical shape of a trapezoid2D.
///
/// It is defined by half lengths in local0 coordinate bounds[0] and bounds[1]
/// at -/+ half length in the local1 coordinate bounds[2]. bounds[3] contains
/// the precomputed value of 1 / (2 * bounds[2]), which avoids
/// excessive floating point divisions.
class trapezoid2D {
 public:
  /// The name for this shape
  static constexpr std::string_view name = "trapezoid2D";

  enum boundaries : unsigned int {
    e_half_length_0 = 0u,
    e_half_length_1 = 1u,
    e_half_length_2 = 2u,
    e_divisor = 3u,  // 1 / (2 * bounds[e_half_length_2])
    e_size = 4u,
  };

  /// Container definition for the shape boundary values
  template <concepts::scalar scalar_t>
  using bounds_type = darray<scalar_t, boundaries::e_size>;

  /// Local coordinate frame for boundary checks
  template <concepts::algebra algebra_t>
  using local_frame_type = cartesian2D<algebra_t>;

  /// Result type of a boundary check
  template <typename bool_t>
  using result_type = detail::boundary_check_result<bool_t>;

  /// Dimension of the local coordinate system
  static constexpr std::size_t dim{2u};

  /// @brief Find the minimum distance to any boundary.
  ///
  /// @note the point is expected to be given in local coordinates by the
  /// caller.
  ///
  /// @param bounds the boundary values for this shape
  /// @param loc_p the point to be checked in the local coordinate system
  ///
  /// @return the minimum distance.
  template <concepts::scalar scalar_t, concepts::point point_t>
  DETRAY_HOST_DEVICE inline scalar_t min_dist_to_boundary(
      const bounds_type<scalar_t> &bounds, const point_t &loc_p) const {
    // Minimum distance between loc point and line between corner points
    const scalar_t d_y = -2.f * bounds[e_half_length_2];
    const scalar_t d_x = bounds[e_half_length_1] - bounds[e_half_length_0];

    const scalar_t denom = math::sqrt(d_y * d_y + d_x * d_x);
    const scalar_t dist_to_line =
        math::fabs(d_y * loc_p[0] - d_x * loc_p[1] -
                   bounds[e_half_length_2] *
                       (bounds[e_half_length_1] + bounds[e_half_length_0])) /
        denom;

    return math::min(dist_to_line, math::fabs(math::fabs(loc_p[1]) -
                                              bounds[e_half_length_2]));
  }

  /// @brief Check boundary values for a local point.
  /// @{
  /// @param bounds the boundary values for this shape
  /// @param trf the surface transform
  /// @param glob_p the point to be checked in the global coordinate system
  /// @param tol dynamic tolerance determined by caller
  ///
  /// @return true if the local point lies within the given boundaries.
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE constexpr result_type<dbool<algebra_t>> check_boundaries(
      const bounds_type<dscalar<algebra_t>> &bounds,
      const dtransform3D<algebra_t> &trf, const dpoint3D<algebra_t> &glob_p,
      const dscalar<algebra_t> tol =
          std::numeric_limits<dscalar<algebra_t>>::epsilon(),
      const dscalar<algebra_t> edge_tol = 0.f) const {
    // Get the full local position
    const dpoint3D<algebra_t> loc_p =
        cartesian3D<algebra_t>::global_to_local(trf, glob_p, {});

    return check_boundaries(bounds, loc_p, tol, edge_tol);
  }

  /// @note the point is expected to be given in local coordinates by the
  /// caller. For the conversion from global cartesian coordinates, the
  /// nested @c shape struct can be used.
  ///
  /// @param bounds the boundary values for this shape
  /// @param loc_p the point to be checked in the local coordinate system
  /// @param tol dynamic tolerance determined by caller
  ///
  /// @return true if the local point lies within the given boundaries.
  template <concepts::scalar scalar_t, concepts::point point_t>
  DETRAY_HOST_DEVICE constexpr auto check_boundaries(
      const bounds_type<scalar_t> &bounds, const point_t &loc_p,
      const scalar_t tol = std::numeric_limits<scalar_t>::epsilon(),
      const scalar_t edge_tol = 0.f) const {
    const scalar_t rel_y =
        (bounds[e_half_length_2] + loc_p[1]) * bounds[e_divisor];
    const scalar_t bound_x{
        bounds[e_half_length_0] +
        rel_y * (bounds[e_half_length_1] - bounds[e_half_length_0])};

    auto inside_mask{(math::fabs(loc_p[0]) <= (bound_x + tol) &&
                      math::fabs(loc_p[1]) <= (bounds[e_half_length_2] + tol))};

    decltype(inside_mask) inside_edge{false};
    if (detail::any_of(edge_tol > 0.f)) {
      const scalar_t full_tol{tol + edge_tol};

      inside_edge =
          (math::fabs(loc_p[0]) <= (bound_x + full_tol) &&
           math::fabs(loc_p[1]) <= (bounds[e_half_length_2] + full_tol));
    }

    return result_type<decltype(inside_mask)>{inside_mask, inside_edge};
  }
  /// @}

  /// @brief Measure of the shape: Area
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the trapezoid area on the plane
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t measure(
      const bounds_type<scalar_t> &bounds) const {
    return area(bounds);
  }

  /// @brief The area of a the shape
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the trapezoid area.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t area(
      const bounds_type<scalar_t> &bounds) const {
    return 2.f * (bounds[e_half_length_0] + bounds[e_half_length_1]) *
           bounds[e_half_length_2];
  }

  /// @brief Merge two trapezoid shapes
  ///
  /// @param bounds the boundary values for this shape
  /// @param o_bounds the boundary values for the other shape
  ///
  /// @returns merged bound values
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr bounds_type<scalar_t> merge(
      const bounds_type<scalar_t> &bounds,
      const bounds_type<scalar_t> &o_bounds) const {
    bounds_type<scalar_t> new_bounds{};

    new_bounds[e_half_length_0] =
        math::max(bounds[e_half_length_0], o_bounds[e_half_length_0]);
    new_bounds[e_half_length_1] =
        math::max(bounds[e_half_length_1], o_bounds[e_half_length_1]);
    new_bounds[e_half_length_2] =
        math::max(bounds[e_half_length_2], o_bounds[e_half_length_2]);

    // Recalculate the memoized divisor
    new_bounds[e_divisor] = 1.f / (2.f * new_bounds[e_half_length_2]);

    return new_bounds;
  }

  /// @brief Lower and upper point for minimal axis aligned bounding box.
  ///
  /// Computes the min and max vertices in a local cartesian frame.
  ///
  /// @param bounds the boundary values for this shape
  /// @param env dynamic envelope around the shape
  ///
  /// @returns and array of coordinates that contains the lower point (first
  /// three values) and the upper point (latter three values) .
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE inline darray<dscalar<algebra_t>, 6> local_min_bounds(
      const bounds_type<dscalar<algebra_t>> &bounds,
      const dscalar<algebra_t> env =
          std::numeric_limits<dscalar<algebra_t>>::epsilon()) const {
    using scalar_t = dscalar<algebra_t>;

    assert(env > 0.f);
    const scalar_t x_bound{(bounds[e_half_length_0] > bounds[e_half_length_1]
                                ? bounds[e_half_length_0]
                                : bounds[e_half_length_1]) +
                           env};
    const scalar_t y_bound{bounds[e_half_length_2] + env};
    return {-x_bound, -y_bound, -env, x_bound, y_bound, env};
  }

  /// @returns the shapes centroid in local cartesian coordinates
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
      const bounds_type<dscalar<algebra_t>> &bounds) const {
    using scalar_t = dscalar<algebra_t>;

    const scalar_t h_2{bounds[e_half_length_2]};
    const scalar_t a_2{bounds[e_half_length_1]};
    const scalar_t b_2{bounds[e_half_length_0]};

    const scalar_t y{2.f * h_2 * (2.f * a_2 + b_2) * 1.f / (3.f * (a_2 + b_2))};

    return {static_cast<scalar_t>(0), y - h_2, static_cast<scalar_t>(0)};
  }

  /// Generate vertices in local cartesian frame
  ///
  /// @param bounds the boundary values for the trapezoid
  /// @param ls is the number of line segments
  ///
  /// @return a generated list of vertices
  template <concepts::algebra algebra_t>
  DETRAY_HOST dvector<dpoint3D<algebra_t>> vertices(
      const bounds_type<dscalar<algebra_t>> &bounds, dindex /*ignored*/) const {
    using point3_t = dpoint3D<algebra_t>;

    // left hand lower corner
    point3_t lh_lc{-bounds[e_half_length_0], -bounds[e_half_length_2],
                   static_cast<dscalar<algebra_t>>(0.f)};
    // right hand lower corner
    point3_t rh_lc{bounds[e_half_length_0], -bounds[e_half_length_2],
                   static_cast<dscalar<algebra_t>>(0.f)};
    // right hand upper corner
    point3_t rh_uc{bounds[e_half_length_1], bounds[e_half_length_2],
                   static_cast<dscalar<algebra_t>>(0.f)};
    // left hand upper corner
    point3_t lh_uc{-bounds[e_half_length_1], bounds[e_half_length_2],
                   static_cast<dscalar<algebra_t>>(0.f)};

    // Return the confining vertices
    return {lh_lc, rh_lc, rh_uc, lh_uc};
  }

  /// @brief Check consistency of boundary values.
  ///
  /// @param bounds the boundary values for this shape
  /// @param os output stream for error messages
  ///
  /// @return true if the bounds are consistent.
  template <concepts::scalar scalar_t>
  DETRAY_HOST constexpr bool check_consistency(
      const bounds_type<scalar_t> &bounds, std::ostream &os) const {
    constexpr auto tol{10.f * std::numeric_limits<scalar_t>::epsilon()};

    if (bounds[e_half_length_0] < tol || bounds[e_half_length_1] < tol) {
      os << "DETRAY ERROR (HOST): Half length in x must be in the range "
            "(0, "
            "numeric_max)"
         << std::endl;
      return false;
    }
    if (bounds[e_half_length_2] < tol) {
      os << "DETRAY ERROR (HOST): Half length in y must be in the range "
            "(0, "
            "numeric_max)"
         << std::endl;
      return false;
    }

    if (const auto div{1.f / (2.f * bounds[e_half_length_2])};
        math::fabs(bounds[e_divisor] - div) > tol) {
      os << "DETRAY ERROR (HOST): Divisor incorrect. Should be: " << div
         << std::endl;
      return false;
    }

    return true;
  }
};

}  // namespace detray
