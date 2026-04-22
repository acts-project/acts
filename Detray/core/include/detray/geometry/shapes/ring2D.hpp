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
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/cartesian3D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/geometry/detail/shape_utils.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string_view>

namespace detray {

/// @brief Geometrical shape of a closed ring.
///
/// It is defined by the two radii bounds[0] and bounds[1],
/// and can be checked with a tolerance in t[0] and t[1].
class ring2D {
 public:
  /// The name for this shape
  static constexpr std::string_view name = "ring2D";

  enum boundaries : unsigned int {
    e_inner_r = 0u,
    e_outer_r = 1u,
    e_size = 2u,
  };

  /// Container definition for the shape boundary values
  template <concepts::scalar scalar_t>
  using bounds_type = darray<scalar_t, boundaries::e_size>;

  /// Local coordinate frame for boundary checks
  template <concepts::algebra algebra_t>
  using local_frame_type = polar2D<algebra_t>;

  /// Result type of a boundary check
  template <typename bool_t>
  using result_type = bool_t;

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
    return math::min(math::fabs(loc_p[0] - bounds[e_inner_r]),
                     math::fabs(bounds[e_outer_r] - loc_p[0]));
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
      const dscalar<algebra_t> /*edge_tol*/ = 0.f) const {
    // Rotate to the local cartesian frame
    const dpoint3D<algebra_t> loc_p =
        cartesian3D<algebra_t>::global_to_local(trf, glob_p, {});

    // Only check the radius
    return check_boundaries(
        bounds,
        dpoint2D<algebra_t>{vector::perp(loc_p),
                            static_cast<dscalar<algebra_t>>(0.f)},
        tol);
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
      const scalar_t /*edge_tol*/ = 0.f) const {
    return ((loc_p[0] + tol) >= bounds[e_inner_r] &&
            loc_p[0] <= (bounds[e_outer_r] + tol));
  }
  /// @}

  /// @brief Measure of the shape: Area
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the ring area on the plane
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t measure(
      const bounds_type<scalar_t> &bounds) const {
    return area(bounds);
  }

  /// @brief The area of a the shape
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the ring area.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t area(
      const bounds_type<scalar_t> &bounds) const {
    return (bounds[e_outer_r] * bounds[e_outer_r] -
            bounds[e_inner_r] * bounds[e_inner_r]) *
           constant<scalar_t>::pi;
  }

  /// @brief Merge two ring/disc shapes
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

    new_bounds[e_inner_r] = math::min(bounds[e_inner_r], o_bounds[e_inner_r]);
    new_bounds[e_outer_r] = math::max(bounds[e_outer_r], o_bounds[e_outer_r]);

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
    const dscalar<algebra_t> r_bound{env + bounds[e_outer_r]};
    return {-r_bound, -r_bound, -env, r_bound, r_bound, env};
  }

  /// @returns the shapes centroid in local cartesian coordinates
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
      const bounds_type<dscalar<algebra_t>> &) const {
    return {0.f, 0.f, 0.f};
  }

  /// Generate vertices in local cartesian frame
  ///
  /// @param bounds the boundary values for the stereo annulus
  /// @param n_seg is the number of line segments
  ///
  /// @return a generated list of vertices
  template <concepts::algebra algebra_t>
  DETRAY_HOST dvector<dpoint3D<algebra_t>> vertices(
      const bounds_type<dscalar<algebra_t>> &, dindex) const {
    throw std::runtime_error("Vertex generation for rings is not implemented");
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
    constexpr auto tol{10.f * std::numeric_limits<scalar_t>::epsilon()};

    if (math::signbit(bounds[e_inner_r]) || bounds[e_outer_r] < tol) {
      os << "DETRAY ERROR (HOST): Radius must be in the range [0, "
            "numeric_max)"
         << std::endl;
      return false;
    }
    if (bounds[e_inner_r] >= bounds[e_outer_r] ||
        math::fabs(bounds[e_inner_r] - bounds[e_outer_r]) < tol) {
      os << "DETRAY ERROR (HOST): Inner radius must be smaller outer "
            "radius.";
      return false;
    }

    return true;
  }
};

}  // namespace detray
