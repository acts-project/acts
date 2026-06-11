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
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/geometry/shapes/cylinder3D.hpp"
#include "detray/navigation/intersection/bounding_box/cuboid_intersector.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/tracks/ray.hpp"

// System include(s)
#include <cassert>
#include <limits>
#include <type_traits>
#include <vector>

namespace detray {

/// An axis aligned bounding box of a given @tparam shape_t
template <typename shape_t, concepts::algebra algebra_t>
class axis_aligned_bounding_volume {
 public:
  /// Define geometric properties
  /// @{
  using scalar_t = dscalar<algebra_t>;
  using shape = shape_t;
  using boundaries = typename shape_t::boundaries;
  template <typename A>
  using local_frame = typename shape_t::template local_frame_type<A>;

  static constexpr std::size_t dim{shape_t::dim};
  /// @}

  /// Default constructor builds an infinite box
  constexpr axis_aligned_bounding_volume() = default;

  /// Constructor from mask boundary values
  template <typename... Args>
  DETRAY_HOST_DEVICE explicit constexpr axis_aligned_bounding_volume(
      std::size_t box_id, Args&&... args)
      : m_mask(box_id, std::forward<Args>(args)...) {}

  /// Construct around an arbitrary surface @param mask
  template <typename mask_t, typename S = shape_t>
    requires std::is_same_v<S, cuboid3D>
  DETRAY_HOST_DEVICE constexpr axis_aligned_bounding_volume(
      const mask_t& mask, std::size_t box_id, const scalar_t envelope)
      : m_mask{mask.local_min_bounds(envelope).values(), box_id} {
    // Make sure the box is actually 'bounding'
    assert(envelope >= std::numeric_limits<scalar_t>::epsilon());
  }

  /// Construct from mask boundary vector
  DETRAY_HOST axis_aligned_bounding_volume(const std::vector<scalar_t>& values,
                                           std::size_t box_id)
      : m_mask(values, box_id) {
    assert(values.size() == shape::boundaries::e_size &&
           " Given number of boundaries does not match mask shape.");
  }

  /// Construct a bounding box around a set of boxes
  /// @note the given bounding volumes need to be defnined in the same
  /// local coordinate system!
  template <typename other_shape_t>
    requires std::is_same_v<
        typename shape::template local_frame_type<algebra_t>,
        typename other_shape_t::template local_frame_type<algebra_t>>
  DETRAY_HOST constexpr axis_aligned_bounding_volume(
      const std::vector<axis_aligned_bounding_volume<other_shape_t, algebra_t>>&
          aabbs,
      std::size_t box_id, const scalar_t env) {
    using loc_point_t = darray<scalar_t, other_shape_t::dim>;

    // Find min/max extent of the local aabb in local coordinates
    constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
    scalar_t min_x{inv};
    scalar_t min_y{inv};
    scalar_t min_z{inv};
    scalar_t max_x{-inv};
    scalar_t max_y{-inv};
    scalar_t max_z{-inv};
    for (const auto& vol : aabbs) {
      const auto min_point = vol.template loc_min<loc_point_t>();
      const auto max_point = vol.template loc_max<loc_point_t>();

      // Check every coordinate of the points
      min_x = min_point[0] < min_x ? min_point[0] : min_x;
      min_y = min_point[1] < min_y ? min_point[1] : min_y;

      max_x = max_point[0] > max_x ? max_point[0] : max_x;
      max_y = max_point[1] > max_y ? max_point[1] : max_y;

      if constexpr (min_point.size() > 2) {
        min_z = min_point[2] < min_z ? min_point[2] : min_z;
        max_z = max_point[2] > max_z ? max_point[2] : max_z;
      }
    }
    m_mask = mask<shape, algebra_t, std::size_t>{
        box_id,      min_x - env, min_y - env, min_z - env,
        max_x + env, max_y + env, max_z + env};
  }

  /// Subscript operator @returns a single box boundary.
  DETRAY_HOST_DEVICE
  constexpr auto operator[](const std::size_t i) const -> scalar_t {
    return m_mask[i];
  }

  /// @returns the bounds of the box, depending on its shape
  DETRAY_HOST_DEVICE
  constexpr auto id() const -> std::size_t { return m_mask.volume_link(); }

  /// @returns the bounds of the box, depending on its shape
  DETRAY_HOST_DEVICE
  constexpr auto bounds() const -> const mask<shape, algebra_t, std::size_t>& {
    return m_mask;
  }

  /// @returns the minimum bounds of the volume in local coordinates
  template <typename point_t>
  DETRAY_HOST_DEVICE constexpr auto loc_min() const -> point_t {
    if constexpr (std::is_same_v<shape, cuboid3D>) {
      return {m_mask[cuboid3D::e_min_x], m_mask[cuboid3D::e_min_y],
              m_mask[cuboid3D::e_min_z]};
    } else if constexpr (std::is_same_v<shape, cylinder3D>) {
      return {-m_mask[cylinder3D::e_max_r], m_mask[cylinder3D::e_min_phi],
              m_mask[cylinder3D::e_min_z]};
    }

    // If the volume shape is not supported, return universal minimum
    assert(false);
    constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
    return point_t{-inv, -inv, -inv};
  }

  /// @returns the maximum bounds of the volume in local coordinates
  template <typename point_t>
  DETRAY_HOST_DEVICE constexpr auto loc_max() const -> point_t {
    if constexpr (std::is_same_v<shape, cuboid3D>) {
      return {m_mask[cuboid3D::e_max_x], m_mask[cuboid3D::e_max_y],
              m_mask[cuboid3D::e_max_z]};
    } else if constexpr (std::is_same_v<shape, cylinder3D>) {
      return {m_mask[cylinder3D::e_max_r], m_mask[cylinder3D::e_max_phi],
              m_mask[cylinder3D::e_max_z]};
    }

    // If the volume shape is not supported, return universal minimum
    // (or compilation error for 2D point)
    assert(false);
    constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
    return point_t{inv, inv, inv};
  }
  /// @returns the minimum bounds of the volume in global cartesian
  /// coordinates
  template <concepts::transform3D transform3_t>
  DETRAY_HOST_DEVICE constexpr auto glob_min(const transform3_t& trf) const ->
      typename transform3_t::point3 {
    using point3_t = typename transform3_t::point3;

    if constexpr (std::is_same_v<shape, cuboid3D>) {
      return trf.point_to_global(loc_min<point3_t>());
    } else if constexpr (std::is_same_v<shape, cylinder3D>) {
      return cylindrical3D<transform3_t>{}.local_to_global(trf, m_mask,
                                                           loc_min<point3_t>());
    }

    // If the volume shape is not supported, return universal minimum
    // (or compilation error for 2D point)
    assert(false);
    constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
    return point3_t{-inv, -inv, -inv};
  }

  /// @returns the maximum bounds of the volume in global cartesian
  /// coordinates
  template <concepts::transform3D transform3_t>
  DETRAY_HOST_DEVICE constexpr auto glob_max(const transform3_t& trf) const ->
      typename transform3_t::point3 {
    using point3_t = typename transform3_t::point3;

    if constexpr (std::is_same_v<shape, cuboid3D>) {
      return trf.point_to_global(loc_max<point3_t>());
    } else if constexpr (std::is_same_v<shape, cylinder3D>) {
      return cylindrical3D<transform3_t>{}.local_to_global(trf, m_mask,
                                                           loc_max<point3_t>());
    }

    // If the volume shape is not supported, return universal minimum
    assert(false);
    constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
    return point3_t{inv, inv, inv};
  }

  /// @returns the geometric center position in global cartesian system
  template <concepts::point3D point3_t>
  DETRAY_HOST_DEVICE constexpr auto center() const -> point3_t {
    const scalar_t center_x{
        0.5f * (m_mask[cuboid3D::e_max_x] + m_mask[cuboid3D::e_min_x])};
    const scalar_t center_y{
        0.5f * (m_mask[cuboid3D::e_max_y] + m_mask[cuboid3D::e_min_y])};
    const scalar_t center_z{
        0.5f * (m_mask[cuboid3D::e_max_z] + m_mask[cuboid3D::e_min_z])};

    return {detail::is_invalid_value(center_x) ? 0.f : center_x,
            detail::is_invalid_value(center_y) ? 0.f : center_y,
            detail::is_invalid_value(center_z) ? 0.f : center_z};
  }

  /// @brief Lower and upper point for minimum axis aligned bounding box of
  /// cuboid shape.
  ///
  /// Computes the min and max vertices in global coordinates from the local
  /// minimum aabb. The global aabb is not necessarily an minimum aabb.
  ///
  /// @param trf affine transformation
  ///
  /// @returns a new, transformed aabb.
  template <concepts::transform3D transform3_t, typename S = shape_t>
    requires std::is_same_v<S, cuboid3D>
  DETRAY_HOST_DEVICE auto transform(const transform3_t& trf) const
      -> axis_aligned_bounding_volume {
    using point3_t = typename transform3_t::point3;

    const scalar_t scalor_x{
        (m_mask[cuboid3D::e_max_x] - m_mask[cuboid3D::e_min_x])};
    const scalar_t scalor_y{
        (m_mask[cuboid3D::e_max_y] - m_mask[cuboid3D::e_min_y])};
    const scalar_t scalor_z{
        (m_mask[cuboid3D::e_max_z] - m_mask[cuboid3D::e_min_z])};

    // Cannot handle 'inv' propagation through the calculation for now
    if (detail::is_invalid_value(scalor_x) ||
        detail::is_invalid_value(scalor_y) ||
        detail::is_invalid_value(scalor_z)) {
      // If the box was infinite to begin with, it stays that way
      assert(detail::is_invalid_value(scalor_x) &&
             detail::is_invalid_value(scalor_y) &&
             detail::is_invalid_value(scalor_z));

      return *this;
    }

    // The new axis vectors, scaled to the aabb dimensions
    // (e.g. max_x - min_x)
    const point3_t new_box_x = scalor_x * trf.x();
    const point3_t new_box_y = scalor_y * trf.y();
    const point3_t new_box_z = scalor_z * trf.z();

    // Transform the old min and max points to the global frame and
    // construct all corner points of the local aabb in global coordinates
    darray<point3_t, 8> glob_c_points;
    glob_c_points[0] = glob_min(trf);
    glob_c_points[1] = glob_max(trf);
    glob_c_points[2] = glob_c_points[0] + new_box_x;
    glob_c_points[3] = glob_c_points[0] + new_box_y;
    glob_c_points[4] = glob_c_points[0] + new_box_z;
    glob_c_points[5] = glob_c_points[1] - new_box_x;
    glob_c_points[6] = glob_c_points[1] - new_box_y;
    glob_c_points[7] = glob_c_points[1] - new_box_z;

    // Find min/max extent of the local aabb in global coordinates
    constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
    scalar_t min_x{inv};
    scalar_t min_y{inv};
    scalar_t min_z{inv};
    scalar_t max_x{-inv};
    scalar_t max_y{-inv};
    scalar_t max_z{-inv};
    for (const point3_t& p : glob_c_points) {
      // Check every coordinate of the point
      min_x = p[0] < min_x ? p[0] : min_x;
      max_x = p[0] > max_x ? p[0] : max_x;
      min_y = p[1] < min_y ? p[1] : min_y;
      max_y = p[1] > max_y ? p[1] : max_y;
      min_z = p[2] < min_z ? p[2] : min_z;
      max_z = p[2] > max_z ? p[2] : max_z;
    }

    // Construct the transformed aabb
    return axis_aligned_bounding_volume{
        m_mask.volume_link(), min_x, min_y, min_z, max_x, max_y, max_z};
  }

  /// Checks whether a point lies inside the box. The point has to be defined
  /// in the coordinate frame that is spanned by the box axes.
  template <concepts::point point_t>
  DETRAY_HOST_DEVICE constexpr bool is_inside(
      const point_t& loc_p,
      const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const {
    return m_mask.is_inside(loc_p, t);
  }

  /// Intersect the box with a ray
  DETRAY_HOST_DEVICE
  constexpr bool intersect(
      const detail::ray<algebra_t>& ray,
      const scalar_t t = std::numeric_limits<scalar_t>::epsilon()) const {
    static_assert(std::is_same_v<shape, cuboid3D>,
                  "aabbs are only implemented in cuboid shape for now");
    return cuboid_intersector{}(ray, m_mask, t);
  }

 private:
  /// Print the bounding volume
  DETRAY_HOST friend std::ostream& operator<<(
      std::ostream& os, const axis_aligned_bounding_volume& aabb) {
    return os << aabb.bounds().to_string();
  }
  /// Keeps the box boundary values and id
  mask<shape, algebra_t, std::size_t> m_mask;
};

}  // namespace detray
