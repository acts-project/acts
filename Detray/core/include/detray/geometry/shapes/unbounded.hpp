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
#include "detray/geometry/detail/shape_utils.hpp"
#include "detray/utils/string_helpers.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string>
#include <string_view>

namespace detray {

/// @brief Wraps any shape, but does not enforce boundaries
template <typename shape_t>
class unbounded {
 public:
  /// Wrapped shape whose boundary checks always accept a point.
  using shape = shape_t;
  /// Boundary index enum inherited from the wrapped shape.
  using boundaries = typename shape::boundaries;

  /// Container definition for the shape boundary values
  template <concepts::scalar scalar_t>
  using bounds_type = darray<scalar_t, boundaries::e_size>;

  /// Convenience member to construct the name
  static constexpr std::string_view name_prefix = "unbounded ";

  /// The name for this shape
  static constexpr utils::string_view_concat2 name{name_prefix, shape::name};

  /// Local coordinate frame for boundary checks
  template <concepts::algebra algebra_t>
  using local_frame_type = typename shape::template local_frame_type<algebra_t>;

  /// Result type of a boundary check
  template <typename bool_t>
  using result_type = detail::boundary_check_result<bool_t>;

  /// Dimension of the local coordinate system
  static constexpr std::size_t dim{shape_t::dim};

  /// @brief Find the minimum distance to any boundary.
  ///
  /// @note the point is expected to be given in local coordinates by the
  /// caller.
  ///
  /// @returns The maximum scalar value because this shape has no active edge.
  template <concepts::scalar scalar_t, concepts::point point_t>
  DETRAY_HOST_DEVICE constexpr scalar_t min_dist_to_boundary(
      const bounds_type<scalar_t>& /*bounds*/, const point_t& /*loc_p*/) const {
    return std::numeric_limits<scalar_t>::max();
  }

  /// @brief Check boundary values for a local point.
  ///
  /// @tparam bounds_t any type of boundary values
  ///
  /// @note the parameters are ignored
  ///
  /// @return always true
  /// @{
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE constexpr result_type<dbool<algebra_t>> check_boundaries(
      const bounds_type<dscalar<algebra_t>>& /*bounds*/,
      const dtransform3D<algebra_t>& /*trf*/,
      const dpoint3D<algebra_t>& /*glob_p*/,
      const dscalar<algebra_t> /*tol*/ = 0.f,
      const dscalar<algebra_t> /*edge_tol*/ = 0.f) const {
    if constexpr (std::same_as<result_type<dbool<algebra_t>>,
                               dbool<algebra_t>>) {
      return true;
    } else {
      return {true, true};
    }
  }

  /// Return an inside result for a point supplied in local coordinates.
  /// @param tol Tolerance argument retained for shape interface compatibility.
  /// @param edge_tol Edge tolerance argument retained for interface compatibility.
  /// @returns An inside result, including the edge flag when the result type has one.
  template <typename bounds_t, concepts::point point_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr auto check_boundaries(
      const bounds_t& /*bounds*/, const point_t& /*loc_p*/,
      const scalar_t tol = 0.f, const scalar_t edge_tol = 0.f) const {
    using bool_t = decltype(tol < edge_tol);

    if constexpr (std::same_as<result_type<bool_t>, bool_t>) {
      return bool_t{true};
    } else {
      return result_type<bool_t>{true, true};
    }
  }
  /// @}

  /// @brief Measure of the shape: Inf
  ///
  /// @param bounds Boundary values accepted for interface compatibility.
  /// @returns Inf.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t measure(
      const bounds_type<scalar_t>& bounds) const {
    if constexpr (dim == 2) {
      return area(bounds);
    } else {
      return volume(bounds);
    }
  }

  /// @brief The area of a the shape
  ///
  /// @returns Inf.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t area(
      const bounds_type<scalar_t>& /*unused*/) const {
    return std::numeric_limits<scalar_t>::max();
  }

  /// @brief The volume of a the shape
  ///
  /// @returns Inf.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t volume(
      const bounds_type<scalar_t>& /*unused*/) const {
    return std::numeric_limits<scalar_t>::max();
  }

  /// @brief Merge two unbounded shapes
  ///
  /// @returns merged bound values
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr bounds_type<scalar_t> merge(
      const bounds_type<scalar_t>& /*bounds*/,
      const bounds_type<scalar_t>& /*o_bounds*/) const {
    return {};
  }

  /// @brief Lower and upper point for minimal axis aligned bounding box.
  ///
  /// Computes the min and max vertices in a local cartesian frame of the
  /// shape it wraps.
  ///
  /// @note The @c check_boundaries method will return 'inside' for points
  /// that are outside the local min bounds! This results in the bounding box
  /// to behave reasonably in a BVH, but still be always intersected
  /// successfully when queried directly.
  ///
  /// @returns and array of coordinates that contains the lower point (first
  /// three values) and the upper point (latter three values).
  /// @param bounds Boundary values forwarded to the wrapped shape.
  /// @param env Envelope forwarded to the wrapped shape.
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE inline darray<dscalar<algebra_t>, 6> local_min_bounds(
      const bounds_type<dscalar<algebra_t>>& bounds,
      const dscalar<algebra_t> env =
          std::numeric_limits<dscalar<algebra_t>>::epsilon()) const {
    return shape{}.template local_min_bounds<algebra_t>(bounds, env);
  }

  /// @param bounds Boundary values forwarded to the wrapped shape.
  /// @returns The wrapped shape's centroid in local Cartesian coordinates.
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE auto centroid(
      const bounds_type<dscalar<algebra_t>>& bounds) const {
    return shape{}.template centroid<algebra_t>(bounds);
  }

  /// Generate vertices in local cartesian frame
  ///
  /// @param bounds Boundary values passed to the wrapped shape.
  /// @param n_seg Number of segments used for curved boundary tessellation.
  ///
  /// @return a generated list of vertices
  template <concepts::algebra algebra_t>
  DETRAY_HOST dvector<dpoint3D<algebra_t>> vertices(
      const bounds_type<dscalar<algebra_t>>& bounds, dindex n_seg) const {
    return shape{}.template vertices<algebra_t>(bounds, n_seg);
  }

  /// @brief Check consistency of boundary values.
  ///
  /// @return true if the bounds are consistent.
  template <concepts::scalar scalar_t>
  DETRAY_HOST constexpr bool check_consistency(
      const bounds_type<scalar_t>& /*bounds*/,
      const std::ostream& /*os*/) const {
    return true;
  }
};

}  // namespace detray
