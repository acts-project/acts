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
#include "detray/tracks/ray.hpp"
#include "detray/utils/invalid_values.hpp"

namespace detray {

/// @brief Checks if a cuboid bonuding box is intersected by a ray.
///
/// This will not compute the intersection point and is meant for aabb
/// intersection instead of navigation.
/// For the implementation see:
/// https://tavianator.com/2022/ray_box_boundary.html
/// and
/// https://github.com/acts-project/acts/blob/85a67758f79f3b01aed4025592ab6a5e5ffbd323/Core/include/Acts/Utilities/BoundingBox.ipp
struct cuboid_intersector {
  /// Operator to find intersections between ray and a cuboid aabb
  ///
  /// @tparam algebra_t is the linear algebra implementation used
  ///
  /// @param ray is the input ray trajectory
  /// @param box is the input mask
  /// @param mask_tolerance is the tolerance for mask edges
  ///
  /// @return true if the cuboid aabb was intersected
  template <concepts::algebra algebra_t, typename mask_t>
  DETRAY_HOST_DEVICE bool operator()(
      const detail::ray<algebra_t> &ray, const mask_t &box,
      const dscalar<algebra_t> /*mask_tolerance*/ = 0.f) const {
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using boundaries = typename mask_t::boundaries;

    const point3_type &ro = ray.pos();
    const vector3_type &rd = ray.dir();
    // @TODO: put vector-vector operator/ in algebra-plugins
    constexpr scalar_type inv{detail::invalid_value<scalar_type>()};
    const vector3_type inv_dir{rd[0] == 0.f ? inv : 1.f / rd[0],
                               rd[1] == 0.f ? inv : 1.f / rd[1],
                               rd[2] == 0.f ? inv : 1.f / rd[2]};

    // This is prob. slow -> @todo refactor masks to hold custom mask values
    const vector3_type min{box[boundaries::e_min_x], box[boundaries::e_min_y],
                           box[boundaries::e_min_z]};
    const vector3_type max{box[boundaries::e_max_x], box[boundaries::e_max_y],
                           box[boundaries::e_max_z]};
    // Find tmin and tmax, which define the segment of the ray that
    // intersects the box
    vector3_type t1 = (min - ro);
    vector3_type t2 = (max - ro);

    /// @TODO: add operator* for two vectors in algebra-plugins
    for (unsigned int i{0u}; i < 3u; ++i) {
      t1[i] *= inv_dir[i];
      t2[i] *= inv_dir[i];
    }

    const bool t_comp = t1[0] > t2[0];
    scalar_type tmin{t_comp ? t2[0] : t1[0]};
    scalar_type tmax{t_comp ? t1[0] : t2[0]};

    for (unsigned int i{0u}; i < 2u; ++i) {
      if (t1[i] > t2[i]) {
        tmin = t2[i] < tmin ? tmin : t2[i];
        tmax = t1[i] > tmax ? tmax : t1[i];
      } else {
        tmin = t1[i] < tmin ? tmin : t1[i];
        tmax = t2[i] > tmax ? tmax : t2[i];
      }
    }

    // Was a valid intersection found ?
    return tmin < tmax;
  }
};

}  // namespace detray
