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

/// @brief Underlying geometry for a single parameter bound mask
///
/// @tparam kCheckIndex is the index of the local point on which the mask is
///         applied
template <unsigned int kCheckIndex = 0u>
class single3D {
 public:
  /// The name for this shape
  static constexpr std::string_view name = "single3D";

  enum boundaries : unsigned int {
    e_lower = 0u,
    e_upper = 1u,
    e_size = 2u,
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
  static constexpr std::size_t dim{1u};

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
    return math::min(math::fabs(loc_p[kCheckIndex] - bounds[e_lower]),
                     math::fabs(bounds[e_upper] - loc_p[kCheckIndex]));
  }

  /// @brief Check boundary values for a local point.
  ///
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
    auto inside_mask{(bounds[e_lower] - tol <= loc_p[kCheckIndex] &&
                      loc_p[kCheckIndex] <= bounds[e_upper] + tol)};

    decltype(inside_mask) inside_edge{false};
    if (detail::any_of(edge_tol > 0.f)) {
      const scalar_t full_tol{tol + edge_tol};

      inside_edge = (bounds[e_lower] - full_tol <= loc_p[kCheckIndex] &&
                     loc_p[kCheckIndex] <= bounds[e_upper] + full_tol);
    }

    return result_type<decltype(inside_mask)>{inside_mask, inside_edge};
  }
  /// @}

  /// @brief Measure of the shape: Range
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the range.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t measure(
      const bounds_type<scalar_t> &bounds) const {
    return area(bounds);
  }

  /// @brief The area of a the shape
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the range.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t area(
      const bounds_type<scalar_t> &bounds) const {
    return math::fabs(bounds[e_upper] - bounds[e_lower]);
  }

  /// @brief Merge two 'single' shapes
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

    new_bounds[e_lower] = math::min(bounds[e_lower], o_bounds[e_lower]);
    new_bounds[e_upper] = math::max(bounds[e_upper], o_bounds[e_upper]);

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
    assert(env > 0.f);
    darray<dscalar<algebra_t>, 6> o_bounds{-env, -env, -env, env, env, env};
    o_bounds[kCheckIndex] += bounds[e_lower];
    o_bounds[3u + kCheckIndex] += bounds[e_upper];
    return o_bounds;
  }

  /// @returns the shapes centroid in local cartesian coordinates
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE auto centroid(
      const bounds_type<dscalar<algebra_t>> &bounds) const {
    using point3_t = dpoint3D<algebra_t>;

    point3_t centr{0.f, 0.f, 0.f};
    centr[kCheckIndex] = 0.5f * (bounds[e_lower] + bounds[e_upper]);

    return centr;
  }

  /// Generate vertices in local cartesian frame
  ///
  /// @param bounds the boundary values for the single value
  /// @param n_seg is the number of line segments
  ///
  /// @return a generated list of vertices
  template <concepts::algebra algebra_t>
  DETRAY_HOST dvector<dpoint3D<algebra_t>> vertices(
      const bounds_type<dscalar<algebra_t>> &, dindex) const {
    throw std::runtime_error(
        "Vertex generation for single value shapes is not implemented");
    return {};
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
    if (bounds[e_upper] < bounds[e_lower]) {
      os << "DETRAY ERROR (HOST): Upper bounds must be smaller than "
            "lower bounds ";
      return false;
    }

    return true;
  }
};

}  // namespace detray
