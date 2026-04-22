// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/cartesian3D.hpp"
#include "detray/geometry/coordinates/line2D.hpp"
#include "detray/geometry/detail/shape_utils.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string_view>
#include <type_traits>

namespace detray {

/// @brief Geometrical shape of a line surface.
///
/// @tparam kSquareCrossSect determines whether the line has a cricular or
///         square cross section. This also changes the local coord. frame.m
///
/// The line can either have a circular or a square cross section. In the first
/// case bounds[0] refers to the radius, while in the second case it is the
/// half length of the square. The second boundary bounds[1] is the half length
/// in z.
template <bool kSquareCrossSect = false>
class line {
 public:
  /// The name for this shape
  static constexpr std::string_view name = "line";

  /// Geometrical cross section of the line
  static constexpr bool square_cross_sect = kSquareCrossSect;

  enum boundaries : unsigned int {
    e_cross_section = 0u,
    e_half_z = 1u,
    e_size = 2u
  };

  /// Container definition for the shape boundary values
  template <concepts::scalar scalar_t>
  using bounds_type = darray<scalar_t, boundaries::e_size>;

  /// Local coordinate frame for boundary checks
  template <concepts::algebra algebra_t>
  using local_frame_type = line2D<algebra_t>;

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
    if constexpr (square_cross_sect) {
      const scalar_t dist_x = math::fabs(
          math::fabs(loc_p[0] * math::cos(loc_p[2])) - bounds[e_cross_section]);
      const scalar_t dist_y = math::fabs(
          math::fabs(loc_p[0] * math::sin(loc_p[2])) - bounds[e_cross_section]);
      const scalar_t dist_z =
          math::fabs(math::fabs(loc_p[1]) - bounds[e_half_z]);

      return math::min(math::min(dist_x, dist_y), dist_z);

    } else {
      return math::min(
          math::fabs(math::fabs(loc_p[0]) - bounds[e_cross_section]),
          math::fabs(math::fabs(loc_p[1]) - bounds[e_half_z]));
    }
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
    const auto loc_3D{cartesian3D<algebra_t>::global_to_local(trf, glob_p, {})};

    if constexpr (square_cross_sect) {
      using scalar_t = dscalar<algebra_t>;

      const scalar_t loc_0{math::fabs(loc_3D[0])};
      const scalar_t loc_1{math::fabs(loc_3D[1])};
      const scalar_t loc_2{math::fabs(loc_3D[2])};

      // Check in local cartesian coordinates instead of line coordinates
      auto inside_mask{(loc_0 <= (bounds[e_cross_section] + tol) &&
                        loc_1 <= (bounds[e_cross_section] + tol) &&
                        loc_2 <= (bounds[e_half_z] + tol))};

      decltype(inside_mask) inside_edge{false};
      if (detail::any_of(edge_tol > 0.f)) {
        const scalar_t full_tol{tol + edge_tol};

        inside_edge = (loc_0 <= (bounds[e_cross_section] + full_tol) &&
                       loc_1 <= (bounds[e_cross_section] + full_tol) &&
                       loc_2 <= (bounds[e_half_z] + full_tol));
      }

      return result_type<dbool<algebra_t>>{inside_mask, inside_edge};
    } else {
      // Check in local 2D line coordinates (sign not needed)
      return check_boundaries(
          bounds, dpoint2D<algebra_t>{vector::perp(loc_3D), loc_3D[2]}, tol,
          edge_tol);
    }
  }

  /// @note the point is expected to be given in local coordinates by the
  /// caller. For the conversion from global cartesian coordinates, the
  /// nested @c shape struct can be used. The point is assumed to be in
  /// the cylinder 2D frame (sgn * r, z, phi) or (sgn * r, z), the latter only
  /// works for the circular cross section line shape.
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
    // For a square cross section (e.g. a cell of drift chamber), we check
    // if (1) x and y of the local cart. point is less than the half cell
    // size and (2) the distance to the point of closest approach on the
    // line from the line center is less than the half line length
    if constexpr (square_cross_sect) {
      const scalar_t loc_0{math::fabs(loc_p[0] * math::cos(loc_p[2]))};
      const scalar_t loc_1{math::fabs(loc_p[0] * math::sin(loc_p[2]))};
      const scalar_t loc_2{math::fabs(loc_p[1])};

      auto inside_mask{(loc_0 <= (bounds[e_cross_section] + tol) &&
                        loc_1 <= (bounds[e_cross_section] + tol) &&
                        loc_2 <= bounds[e_half_z] + tol)};

      decltype(inside_mask) inside_edge{false};
      if (detail::any_of(edge_tol > 0.f)) {
        const scalar_t full_tol{tol + edge_tol};

        inside_edge = (loc_0 <= (bounds[e_cross_section] + full_tol) &&
                       loc_1 <= (bounds[e_cross_section] + full_tol) &&
                       loc_2 <= bounds[e_half_z] + full_tol);
      }

      return result_type<decltype(inside_mask)>{inside_mask, inside_edge};
    }
    // For a circular cross section (e.g. straw tube), we check if (1) the
    // radial distance is within the scope and (2) the distance to the point
    // of closest approach on the line from the line center is less than the
    // line half length
    else {
      const scalar_t loc_0{math::fabs(loc_p[0])};
      const scalar_t loc_1{math::fabs(loc_p[1])};

      auto inside_mask{(loc_0 <= (bounds[e_cross_section] + tol) &&
                        loc_1 <= (bounds[e_half_z] + tol))};

      decltype(inside_mask) inside_edge{false};
      if (detail::any_of(edge_tol > 0.f)) {
        const scalar_t full_tol{tol + edge_tol};

        inside_edge = (loc_0 <= (bounds[e_cross_section] + full_tol) &&
                       loc_1 <= (bounds[e_half_z] + full_tol));
      }

      return result_type<decltype(inside_mask)>{inside_mask, inside_edge};
    }
  }
  /// @}

  /// @brief Measure of the shape: Volume
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the line volume.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t measure(
      const bounds_type<scalar_t> &bounds) const {
    if constexpr (square_cross_sect) {
      return 8.f * bounds[e_half_z] * bounds[e_cross_section] *
             bounds[e_cross_section];
    } else {
      return constant<scalar_t>::pi * 2.f * bounds[e_half_z] *
             bounds[e_cross_section] * bounds[e_cross_section];
    }
  }

  /// @brief The area of a the shape
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the stereo annulus area.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t area(
      const bounds_type<scalar_t> &bounds) const {
    if constexpr (square_cross_sect) {
      return 16.f * bounds[e_half_z] * bounds[e_cross_section];
    } else {
      return 4.f * constant<scalar_t>::pi * bounds[e_cross_section] *
             bounds[e_half_z];
    }
  }

  /// @brief Merge two line shapes
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

    new_bounds[e_cross_section] =
        math::max(bounds[e_cross_section], o_bounds[e_cross_section]);
    new_bounds[e_half_z] = math::max(bounds[e_half_z], o_bounds[e_half_z]);

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
    const scalar_t xy_bound{bounds[e_cross_section] + env};
    const scalar_t z_bound{bounds[e_half_z] + env};

    return {-xy_bound, -xy_bound, -z_bound, xy_bound, xy_bound, z_bound};
  }

  /// @returns the shapes centroid in local cartesian coordinates
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
      const bounds_type<dscalar<algebra_t>> &) const {
    return {0.f, 0.f, 0.f};
  }

  /// Generate vertices in local cartesian frame
  ///
  /// @param bounds the boundary values for the line
  /// @param n_seg is the number of line segments
  ///
  /// @return a generated list of vertices
  template <concepts::algebra algebra_t>
  DETRAY_HOST dvector<dpoint3D<algebra_t>> vertices(
      const bounds_type<dscalar<algebra_t>> &bounds, dindex /*ignored*/) const {
    using point3_t = dpoint3D<algebra_t>;

    point3_t lc = {0.f, 0.f, -bounds[e_half_z]};
    point3_t rc = {0.f, 0.f, bounds[e_half_z]};

    return {lc, rc};
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

    if (bounds[e_cross_section] < tol) {
      os << "DETRAY ERROR (HOST): Radius/sides must be in the range (0, "
            "numeric_max)"
         << std::endl;
      return false;
    }
    if (bounds[e_half_z] < tol) {
      os << "DETRAY ERROR (HOST): Half length z must be in the range (0, "
            "numeric_max)"
         << std::endl;
      return false;
    }

    return true;
  }
};

// Radial crossection, boundary check in polar coordinates
using line_circular = line<false>;
// Square crossection, boundary check in cartesian coordinates
using line_square = line<true>;

}  // namespace detray
