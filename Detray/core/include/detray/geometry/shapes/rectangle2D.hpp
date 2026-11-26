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
#include "detray/geometry/detail/shape_utils.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string_view>

namespace detray {

/// @brief Geometrical shape of an axis-aligned rectangle2D.
///
/// The rectangle is defined by the absolute positions of its lower left
/// corner (bounds[e_min_x], bounds[e_min_y]) and its upper right corner
/// (bounds[e_max_x], bounds[e_max_y]) in the local cartesian frame of the
/// surface. The rectangle does not have to be centered at the origin of the
/// local frame. Multiple rectangles can therefore share a single surface
/// (and transform).
class rectangle2D {
 public:
  /// The name for this shape
  static constexpr std::string_view name = "rectangle2D";

  /// Indices describing the entries in the shape boundary array.
  /// @note The order is the same as in the ACTS @c RectangleBounds
  enum boundaries : unsigned int {
    e_min_x = 0u,
    e_min_y = 1u,
    e_max_x = 2u,
    e_max_y = 3u,
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
    const scalar_t dist_x{math::min(math::fabs(loc_p[0] - bounds[e_min_x]),
                                    math::fabs(loc_p[0] - bounds[e_max_x]))};
    const scalar_t dist_y{math::min(math::fabs(loc_p[1] - bounds[e_min_y]),
                                    math::fabs(loc_p[1] - bounds[e_max_y]))};

    return math::min(dist_x, dist_y);
  }

  /// @brief Check boundary values for a local point.
  /// @{
  /// @param bounds the boundary values for this shape
  /// @param trf the surface transform
  /// @param glob_p the point to be checked in the global coordinate system
  /// @param tol dynamic tolerance determined by caller
  /// @param edge_tol Additional tolerance used to classify points near an edge.
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
    const dpoint2D<algebra_t> loc_p =
        local_frame_type<algebra_t>::global_to_local(trf, glob_p, {});

    return check_boundaries(bounds, loc_p, tol, edge_tol);
  }

  /// @note the point is expected to be given in local coordinates by the
  /// caller. For the conversion from global cartesian coordinates, the
  /// nested @c shape struct can be used.
  ///
  /// @param bounds the boundary values for this shape
  /// @param loc_p the point to be checked in the local coordinate system
  /// @param tol dynamic tolerance determined by caller
  /// @param edge_tol Additional tolerance used to classify points near an edge.
  ///
  /// @return true if the local point lies within the given boundaries.
  template <concepts::scalar scalar_t, concepts::point point_t>
  DETRAY_HOST_DEVICE constexpr auto check_boundaries(
      const bounds_type<scalar_t> &bounds, const point_t &loc_p,
      const scalar_t tol = std::numeric_limits<scalar_t>::epsilon(),
      const scalar_t edge_tol = 0.f) const {
    const scalar_t loc_0{loc_p[0]};
    const scalar_t loc_1{loc_p[1]};

    auto inside_mask{
        ((bounds[e_min_x] - tol) <= loc_0 && loc_0 <= (bounds[e_max_x] + tol) &&
         (bounds[e_min_y] - tol) <= loc_1 && loc_1 <= (bounds[e_max_y] + tol))};

    decltype(inside_mask) inside_edge{false};
    if (detail::any_of(edge_tol > 0.f)) {
      const scalar_t full_tol{tol + edge_tol};

      inside_edge = ((bounds[e_min_x] - full_tol) <= loc_0 &&
                     loc_0 <= (bounds[e_max_x] + full_tol) &&
                     (bounds[e_min_y] - full_tol) <= loc_1 &&
                     loc_1 <= (bounds[e_max_y] + full_tol));
    }

    return result_type<decltype(inside_mask)>{inside_mask, inside_edge};
  }
  /// @}

  /// @brief Measure of the shape: Area
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the rectangle area on the plane
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t measure(
      const bounds_type<scalar_t> &bounds) const {
    return area(bounds);
  }

  /// @brief The area of a the shape
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the rectangle area.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t area(
      const bounds_type<scalar_t> &bounds) const {
    return (bounds[e_max_x] - bounds[e_min_x]) *
           (bounds[e_max_y] - bounds[e_min_y]);
  }

  /// @brief Merge two rectangle shapes
  ///
  /// The result is the smallest axis aligned rectangle that contains both
  /// rectangles.
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

    new_bounds[e_min_x] = math::min(bounds[e_min_x], o_bounds[e_min_x]);
    new_bounds[e_min_y] = math::min(bounds[e_min_y], o_bounds[e_min_y]);
    new_bounds[e_max_x] = math::max(bounds[e_max_x], o_bounds[e_max_x]);
    new_bounds[e_max_y] = math::max(bounds[e_max_y], o_bounds[e_max_y]);

    return new_bounds;
  }

  /// @brief Lower and upper point for minimal axis aligned bounding box.
  ///
  /// Computes the min and max vertices in a local cartesian frame.
  ///
  /// @param bounds the boundary values for this shape
  /// @param env dynamic envelope around the shape
  ///
  /// @returns an array of coordinates that contains the lower point (first
  /// three values) and the upper point (latter three values) .
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE inline darray<dscalar<algebra_t>, 6> local_min_bounds(
      const bounds_type<dscalar<algebra_t>> &bounds,
      const dscalar<algebra_t> env =
          std::numeric_limits<dscalar<algebra_t>>::epsilon()) const {
    assert(env > 0.f);
    return {bounds[e_min_x] - env, bounds[e_min_y] - env, -env,
            bounds[e_max_x] + env, bounds[e_max_y] + env, env};
  }

  /// @param bounds the boundary values for this shape
  ///
  /// @returns the shapes centroid in local cartesian coordinates
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
      const bounds_type<dscalar<algebra_t>> &bounds) const {
    using scalar_t = dscalar<algebra_t>;

    return {0.5f * (bounds[e_min_x] + bounds[e_max_x]),
            0.5f * (bounds[e_min_y] + bounds[e_max_y]),
            static_cast<scalar_t>(0.f)};
  }

  /// Generate vertices in local cartesian frame
  ///
  /// @param bounds Corner positions of the rectangle.
  /// The second argument is unused because a rectangle always has four
  /// corners.
  ///
  /// @return a generated list of vertices
  template <concepts::algebra algebra_t>
  DETRAY_HOST dvector<dpoint3D<algebra_t>> vertices(
      const bounds_type<dscalar<algebra_t>> &bounds, dindex /*ignored*/) const {
    using scalar_t = dscalar<algebra_t>;
    using point3_t = dpoint3D<algebra_t>;

    constexpr auto z{static_cast<scalar_t>(0.f)};

    // left hand lower corner
    point3_t lh_lc{bounds[e_min_x], bounds[e_min_y], z};
    // right hand lower corner
    point3_t rh_lc{bounds[e_max_x], bounds[e_min_y], z};
    // right hand upper corner
    point3_t rh_uc{bounds[e_max_x], bounds[e_max_y], z};
    // left hand upper corner
    point3_t lh_uc{bounds[e_min_x], bounds[e_max_y], z};

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
    if (constexpr auto tol{10.f * std::numeric_limits<scalar_t>::epsilon()};
        (bounds[e_max_x] - bounds[e_min_x]) < tol ||
        (bounds[e_max_y] - bounds[e_min_y]) < tol) {
      os << "DETRAY ERROR (HOST): The upper corner must be larger than the "
            "lower corner in x and y"
         << std::endl;
      return false;
    }

    return true;
  }
};

}  // namespace detray
