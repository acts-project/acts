// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/material/material.hpp"
#include "detray/navigation/accelerators/search_window.hpp"
#include "detray/navigation/concepts.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/quadratic_equation.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <algorithm>
#include <vector>

namespace detray::detail {

/// A functor to retrieve the volume material parameters at a given position
struct get_material_params {
  template <typename mat_group_t, typename index_t, concepts::point point_t>
  DETRAY_HOST_DEVICE inline auto operator()(const mat_group_t &mat_group,
                                            const index_t &idx,
                                            const point_t &loc_p) const {
    using material_t = typename mat_group_t::value_type;

    if constexpr (concepts::volume_material<material_t>) {
      if constexpr (concepts::homogeneous_material<material_t>) {
        // Homogeneous volume material
        return &(detail::material_accessor::get(mat_group, idx, loc_p));
      } else {
        // Volume material maps
        return &(detail::material_accessor::get(mat_group, idx, loc_p)
                     .get_material());
      }
    } else {
      using scalar_t = typename material_t::scalar_type;
      // Cannot be reached for volumes
      return static_cast<const material<scalar_t> *>(nullptr);
    }
  }
};

/// A functor to access all surfaces registered in an acceleration structure
/// @tparam functor_t functor that performs a task per surface
template <typename functor_t>
struct apply_to_surfaces {
  /// Call operator that forwards the functor to all contained surfaces
  template <concepts::accelerator_collection accel_coll_t,
            typename accel_index_t, typename... Args>
  DETRAY_HOST_DEVICE inline void operator()(const accel_coll_t &coll,
                                            const accel_index_t index,
                                            Args &&...args) const {
    using accel_type = typename accel_coll_t::value_type;

    if constexpr (concepts::surface_accelerator<accel_type>) {
      // Run over the surfaces in a single acceleration data structure
      for (const auto &sf : coll[index].all()) {
        functor_t{}(sf, std::forward<Args>(args)...);
      }
    }
  }
};

/// A functor to find surfaces in the neighborhood of a track position
/// @tparam functor_t functor that performs a task per surface in the
///                   neighbourhood (e.g. intersection)
template <typename functor_t>
struct apply_to_neighbourhood {
  /// Call operator that forwards the neighborhood search call in a volume
  /// to a surface finder data structure
  template <concepts::accelerator_collection accel_coll_t,
            typename accel_index_t, concepts::detector detector_t,
            typename track_t, typename window_size_t, typename... Args>
  DETRAY_HOST_DEVICE inline void operator()(
      const accel_coll_t &coll, const accel_index_t index,
      const detector_t &det, const typename detector_t::volume_type &volume,
      const track_t &track, const search_window<window_size_t, 2> &win_size,
      const typename detector_t::geometry_context &ctx, Args &&...args) const {
    using accel_type = typename accel_coll_t::value_type;

    decltype(auto) accel = coll[index];

    if constexpr (concepts::surface_accelerator<accel_type>) {
      // Run over the surfaces in a single acceleration data structure
      for (const auto &sf : accel.search(det, volume, track, win_size, ctx)) {
        assert(sf.volume() == volume.index());
        functor_t{}(sf, std::forward<Args>(args)...);
      }
    } else if constexpr (concepts::volume_accelerator<accel_type>) {
      // Run over the daughter volumes in a single acceleration data
      // structure
      for (const dindex daughter_idx :
           accel.search(det, volume, track, win_size, ctx)) {
        functor_t{}(daughter_idx, std::forward<Args>(args)...);
      }
    }
  }
};

/// A functor to access the daughter volumes of a volume
///
/// Calls the functor once for every daughter volume, also if the
/// acceleration structure contains the daughter volume more than once (e.g.
/// in several bins of a volume grid). Invalid entries (e.g. empty grid bins)
/// are skipped. Runs on the host only.
template <typename functor_t>
struct apply_to_volumes {
  /// Call operator that forwards the functor call to all contained daughter
  /// volumes
  ///
  /// @param coll the collection of volume acceleration structures
  /// @param index the index of the acceleration structure of the volume
  /// @param args additional arguments for the functor
  template <concepts::accelerator_collection accel_coll_t,
            typename accel_index_t, typename... Args>
  DETRAY_HOST inline void operator()(const accel_coll_t &coll,
                                     const accel_index_t index,
                                     Args &&...args) const {
    using accel_type = typename accel_coll_t::value_type;

    if constexpr (concepts::volume_accelerator<accel_type>) {
      // Unique and valid daughter volume indices
      std::vector<dindex> daughters{};
      for (const dindex daughter_idx : coll[index].all()) {
        if (daughter_idx != dindex_invalid) {
          daughters.push_back(daughter_idx);
        }
      }
      std::ranges::sort(daughters);
      const auto [first, last] = std::ranges::unique(daughters);
      daughters.erase(first, last);

      for (const dindex daughter_idx : daughters) {
        functor_t{}(daughter_idx, args...);
      }
    }
  }
};

/// @brief Collect the bounds of a cylinder volume from its portals
///
/// Is called through the mask store visitor for every portal of a cylinder
/// volume: Discs give the radial bounds, cylinder mantles give the bounds in
/// z. Other portal shapes (e.g. planes for phi sectors) are not supported.
struct cylinder_volume_bounds {
  /// Call operator for the mask store visitor
  ///
  /// @tparam mask_group_t the type of the mask collection
  /// @tparam mask_range_t the type of the mask index (range) of the portal
  /// @tparam scalar_t the scalar type
  ///
  /// @param mask_group the collection of masks of the same type
  /// @param mask_range the index (range) of the masks of the portal
  /// @param z_offset position of the portal on the volume z-axis
  /// @param z_sign +1 if the portal z-axis points along the volume z-axis,
  ///               -1 otherwise
  /// @param bounds r_min, r_max, z_min, z_max of the volume (in/out)
  /// @param is_supported set to false for unsupported portal shapes (out)
  template <typename mask_group_t, typename mask_range_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline void operator()(const mask_group_t &mask_group,
                                            const mask_range_t &mask_range,
                                            const scalar_t z_offset,
                                            const scalar_t z_sign,
                                            darray<scalar_t, 4> &bounds,
                                            bool &is_supported) const {
    using shape_t = typename mask_group_t::value_type::shape;

    for (const auto &mask : detray::ranges::subrange(mask_group, mask_range)) {
      if constexpr (std::same_as<shape_t, ring2D>) {
        bounds[0] = math::min(bounds[0], scalar_t{mask[ring2D::e_inner_r]});
        bounds[1] = math::max(bounds[1], scalar_t{mask[ring2D::e_outer_r]});
      } else if constexpr (std::same_as<shape_t, concentric_cylinder2D> ||
                           std::same_as<shape_t, cylinder2D>) {
        const scalar_t z0{z_offset + z_sign * mask[shape_t::e_lower_z]};
        const scalar_t z1{z_offset + z_sign * mask[shape_t::e_upper_z]};
        bounds[2] = math::min(bounds[2], math::min(z0, z1));
        bounds[3] = math::max(bounds[3], math::max(z0, z1));
      } else {
        is_supported = false;
      }
    }
  }
};

/// @brief Count the crossings of a ray with a portal surface
///
/// Is used to decide whether a point lies inside a volume: A ray from the
/// point crosses the portals of the volume an odd number of times, if the
/// point is inside the volume. Works for planar and cylindrical portals.
///
/// The functor flags the test as ambiguous, if the ray passes close to a mask
/// edge (e.g. a seam between two portal surfaces), if the ray touches a
/// cylinder or if the ray runs inside the surface. The caller should then
/// repeat the test with a different ray direction.
struct count_portal_crossings {
  /// Call operator for the mask store visitor
  ///
  /// @param mask_group the collection of masks of the same type
  /// @param mask_range the index (range) of the masks of the portal
  /// @param ray the test ray
  /// @param trf the placement transform of the portal
  /// @param tol the tolerance for edges and for points on the portal
  /// @param n_crossings the number of crossings (output)
  /// @param is_ambiguous set to true if the test is ambiguous (output)
  /// @param is_on_portal set to true if the ray origin is on the portal
  ///                     (output)
  template <typename mask_group_t, typename mask_range_t,
            concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE inline void operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      const detail::ray<algebra_t> &ray, const dtransform3D<algebra_t> &trf,
      const dscalar<algebra_t> tol, unsigned int &n_crossings,
      bool &is_ambiguous, bool &is_on_portal) const {
    using mask_t = typename mask_group_t::value_type;
    using frame_t = typename mask_t::local_frame;
    using scalar_t = dscalar<algebra_t>;
    using point3_t = dpoint3D<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    constexpr scalar_t eps{1e-6f};

    const point3_t &ro = ray.pos();
    const vector3_t &rd = ray.dir();

    // Path lengths to the underlying geometry (at most two solutions)
    darray<scalar_t, 2> paths{};
    unsigned int n_paths{0u};

    if constexpr (concepts::planar_frame<frame_t>) {
      const vector3_t &normal = trf.z();
      const scalar_t denom{vector::dot(normal, rd)};
      const scalar_t num{vector::dot(normal, trf.translation() - ro)};
      if (math::fabs(denom) < eps) {
        // Ray runs inside the plane
        if (math::fabs(num) < tol) {
          is_ambiguous = true;
        }
        return;
      }
      paths[0] = num / denom;
      n_paths = 1u;
    } else if constexpr (std::same_as<frame_t, cylindrical2D<algebra_t>> ||
                         std::same_as<frame_t,
                                      concentric_cylindrical2D<algebra_t>>) {
      // Use the first mask for the radius (all masks share the cylinder)
      dindex mask_idx{0u};
      if constexpr (concepts::interval<mask_range_t>) {
        mask_idx = mask_range.lower();
      } else {
        mask_idx = mask_range;
      }
      const scalar_t r{mask_group[mask_idx][mask_t::shape::e_r]};

      const vector3_t &sz = trf.z();
      const vector3_t pc_cross_sz = vector::cross(ro - trf.translation(), sz);
      const vector3_t rd_cross_sz = vector::cross(rd, sz);
      const scalar_t a{vector::dot(rd_cross_sz, rd_cross_sz)};
      const scalar_t b{2.f * vector::dot(rd_cross_sz, pc_cross_sz)};
      const scalar_t c{vector::dot(pc_cross_sz, pc_cross_sz) - (r * r)};

      // Ray parallel to the cylinder axis
      if (a < eps) {
        if (math::fabs(vector::norm(pc_cross_sz) - r) < tol) {
          is_ambiguous = true;
        }
        return;
      }

      const detail::quadratic_equation<scalar_t> qe{a, b, c};
      if (qe.solutions() == 1) {
        // The ray touches the cylinder
        is_ambiguous = true;
        return;
      }
      if (qe.solutions() == 2) {
        paths[0] = qe.smaller();
        paths[1] = qe.larger();
        n_paths = 2u;
      }
    } else {
      // Other surface shapes cannot form volume boundaries
      return;
    }

    for (unsigned int i = 0u; i < n_paths; ++i) {
      const scalar_t path{paths[i]};
      if (path < -tol) {
        continue;
      }
      const point3_t glob = ro + path * rd;

      // Check all masks of the portal
      for (const auto &mask :
           detray::ranges::subrange(mask_group, mask_range)) {
        const auto loc = mask.to_local_frame(trf, glob, rd);

        // Hit within tolerance: Check the distance to the mask edges
        if (!mask.is_inside(loc, tol)) {
          continue;
        }
        if constexpr (requires(const typename mask_t::shape &sh) {
                        sh.min_dist_to_boundary(mask.values(), loc);
                      }) {
          const scalar_t dist{
              mask.get_shape().min_dist_to_boundary(mask.values(), loc)};
          if (dist < tol) {
            // Close to an edge: The hit might belong to a neighbour portal
            is_ambiguous = true;
          }
        }

        if (path < tol) {
          // The ray origin lies on the portal
          is_on_portal = true;
        } else {
          ++n_crossings;
        }
        break;
      }
    }
  }
};

}  // namespace detray::detail
