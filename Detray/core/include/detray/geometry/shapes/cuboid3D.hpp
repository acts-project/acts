// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/coordinates/cartesian3D.hpp"
#include "detray/geometry/detail/shape_utils.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string_view>

namespace detray {

/// @brief Geometrical shape of a full 3D cuboid.
///
/// It is defined by 3 min-max length pairs and checks whether a point is
/// somewhere inside the cuboid. This type is mainly used for aabb description.
class cuboid3D {
 public:
  /// The name for this shape
  static constexpr std::string_view name = "cuboid3D";

  enum boundaries : unsigned int {
    e_min_x = 0u,
    e_min_y = 1u,
    e_min_z = 2u,
    e_max_x = 3u,
    e_max_y = 4u,
    e_max_z = 5u,
    e_size = 6u,
  };

  /// Container definition for the shape boundary values
  template <concepts::scalar scalar_t>
  using bounds_type = darray<scalar_t, boundaries::e_size>;

  /// Local coordinate frame for boundary checks
  template <concepts::algebra algebra_t>
  using local_frame_type = cartesian3D<algebra_t>;

  /// Result type of a boundary check
  template <typename bool_t>
  using result_type = bool_t;

  /// Dimension of the local coordinate system
  static constexpr std::size_t dim{3u};

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
    const scalar_t min_x_dist =
        math::min(math::fabs(loc_p[0] - bounds[e_min_x]),
                  math::fabs(bounds[e_max_x] - loc_p[0]));
    const scalar_t min_y_dist =
        math::min(math::fabs(loc_p[1] - bounds[e_min_y]),
                  math::fabs(bounds[e_max_y] - loc_p[1]));
    const scalar_t min_z_dist =
        math::min(math::fabs(loc_p[2] - bounds[e_min_z]),
                  math::fabs(bounds[e_max_z] - loc_p[2]));

    // Use the chord for the phi distance
    return math::min(math::min(min_x_dist, min_y_dist), min_z_dist);
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
      const dscalar<algebra_t> /*edge_tol*/ = 0.f) const {
    // Get the full local position
    const dpoint3D<algebra_t> loc_p =
        local_frame_type<algebra_t>::global_to_local(trf, glob_p, {});

    return check_boundaries(bounds, loc_p, tol);
  }

  /// @note the point is expected to be given in local coordinates by the
  /// caller. For the conversion from global cartesian coordinates, the
  /// nested @c shape struct can be used. The point is assumed to be in
  /// the a cartesian 3D frame.
  ///
  /// @param bounds the boundary values for this shape
  /// @param loc_p the point to be checked in the local coordinate system
  /// @param tol dynamic tolerance determined by caller
  ///
  /// @return true if the local point lies within the given boundaries.
  template <concepts::scalar scalar_t, concepts::point3D point_t>
  DETRAY_HOST_DEVICE constexpr auto check_boundaries(
      const bounds_type<scalar_t> &bounds, const point_t &loc_p,
      const scalar_t tol = std::numeric_limits<scalar_t>::epsilon(),
      const scalar_t /*edge_tol*/ = 0.f) const {
    return ((bounds[e_min_x] - tol) <= loc_p[0] &&
            (bounds[e_min_y] - tol) <= loc_p[1] &&
            (bounds[e_min_x] - tol) <= loc_p[2] &&
            loc_p[0] <= (bounds[e_max_x] + tol) &&
            loc_p[1] <= (bounds[e_max_y] + tol) &&
            loc_p[2] <= (bounds[e_max_z] + tol));
  }
  /// @}

  /// @brief Measure of the shape: Volume
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the cuboid volume as part of global space.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t measure(
      const bounds_type<scalar_t> &bounds) const {
    return volume(bounds);
  }

  /// @brief The volume of a the shape
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the cuboid volume.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t volume(
      const bounds_type<scalar_t> &bounds) const {
    return (bounds[e_max_x] - bounds[e_min_x]) *
           (bounds[e_max_y] - bounds[e_min_y]) *
           (bounds[e_max_z] - bounds[e_min_z]);
  }

  /// @brief Merge two 3D cuboid shapes
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
    new_bounds[e_min_z] = math::min(bounds[e_min_z], o_bounds[e_min_z]);
    new_bounds[e_max_x] = math::max(bounds[e_max_x], o_bounds[e_max_x]);
    new_bounds[e_max_y] = math::max(bounds[e_max_y], o_bounds[e_max_y]);
    new_bounds[e_max_z] = math::max(bounds[e_max_z], o_bounds[e_max_z]);

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
    bounds_type<dscalar<algebra_t>> o_bounds{bounds};
    for (unsigned int i{0}; i < 3u; ++i) {
      o_bounds[i] -= env;
      o_bounds[i + 3u] += env;
    }
    return o_bounds;
  }

  /// @returns the shapes centroid in local cartesian coordinates
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
      const bounds_type<dscalar<algebra_t>> &bounds) const {
    return 0.5f * dpoint3D<algebra_t>{bounds[e_min_x] + bounds[e_max_x],
                                      bounds[e_min_y] + bounds[e_max_y],
                                      bounds[e_min_z] + bounds[e_max_z]};
  }

  /// Generate vertices in local cartesian frame
  ///
  /// @param bounds the boundary values for the cuboid
  /// @param n_seg is the number of line segments
  ///
  /// @return a generated list of vertices
  template <concepts::algebra algebra_t>
  DETRAY_HOST dvector<dpoint3D<algebra_t>> vertices(
      const bounds_type<dscalar<algebra_t>> &, dindex) const {
    throw std::runtime_error(
        "Vertex generation for cuboids is not implemented");
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

    if (bounds[e_min_x] >= bounds[e_max_x] ||
        math::fabs(bounds[e_min_x] - bounds[e_max_x]) < tol) {
      os << "DETRAY ERROR (HOST): Min x must be smaller than max x.";
      return false;
    }
    if (bounds[e_min_y] >= bounds[e_max_y] ||
        math::fabs(bounds[e_min_y] - bounds[e_max_y]) < tol) {
      os << "DETRAY ERROR (HOST): Min y must be smaller than max y.";
      return false;
    }
    if (bounds[e_min_z] >= bounds[e_max_z] ||
        math::fabs(bounds[e_min_z] - bounds[e_max_z]) < tol) {
      os << "DETRAY ERROR (HOST): Min z must be smaller than max z.";
      return false;
    }

    return true;
  }
};

}  // namespace detray
