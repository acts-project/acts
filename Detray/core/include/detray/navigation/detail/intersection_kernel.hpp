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
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::detail {

/// The maximum number of solutions that any intersector can produce
inline constexpr std::uint8_t max_n_solutions{2u};

/// A functor to find all valid intersections between the trajectory and surface
///
/// The intersections are returned to the caller instead of being written into
/// a container directly, so that the container code is instantiated once,
/// rather than once per mask type in the mask store visit.
template <template <typename, typename, bool> class intersector_t,
          typename intersection_t>
struct intersection_initialize_impl {
  /// Operator function to initialize intersections
  ///
  /// @tparam mask_group_t is the input mask group type found by variadic
  /// unrolling
  /// @tparam traj_t is the input trajectory type (e.g. ray or helix)
  /// @tparam surface_t is the input surface type
  /// @tparam transform_container_t is the input transform store type
  ///
  /// @param mask_group is the input mask group
  /// @param traj is the input trajectory
  /// @param surface is the input surface
  /// @param contextual_transforms is the input transform container
  /// @param mask_tolerance is the tolerance for mask size
  /// @param overstep_tol negative cutoff for the path
  ///
  /// @return one slot per solution of the intersector. A slot that holds no
  ///         valid intersection is left default constructed, i.e. 'outside'.
  template <typename mask_group_t, typename mask_range_t, typename traj_t,
            typename surface_t, typename transform_container_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline darray<intersection_t, max_n_solutions> operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      const traj_t &traj, const surface_t &sf_desc,
      const transform_container_t &contextual_transforms,
      const typename transform_container_t::context_type &ctx,
      const intersection::config &cfg,
      const scalar_t external_mask_tolerance = 0.f) const {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;
    using algebra_t = typename mask_t::algebra_type;

    darray<intersection_t, max_n_solutions> found{};

    // Find the point of intersection with the underlying geometry
    const auto &ctf = contextual_transforms.at(sf_desc.transform(), ctx);

    constexpr intersector_t<shape_t, algebra_t, intersection_t::contains_pos()>
        intersector{};

    constexpr std::uint8_t n_sol{decltype(intersector)::n_solutions};
    static_assert(n_sol <= max_n_solutions);

    typename decltype(intersector)::result_type result{};

    if constexpr (concepts::cylindrical<mask_t>) {
      dindex mask_idx{detail::invalid_value<dindex>()};
      if constexpr (concepts::interval<mask_range_t>) {
        mask_idx = mask_range.lower();
      } else {
        mask_idx = mask_range;
      }
      assert(mask_idx < mask_group.size());

      result = intersector.point_of_intersection(
          traj, ctf, mask_group[mask_idx], cfg.overstep_tolerance);
    } else {
      result =
          intersector.point_of_intersection(traj, ctf, cfg.overstep_tolerance);
    }

    // Check if any valid solutions were found
    if constexpr (n_sol > 1) {
      bool found_any{false};
      for (const auto &ip : result) {
        if (ip.is_valid()) {
          found_any = true;
        }
      }
      if (!found_any) [[unlikely]] {
        return found;
      }
    } else {
      if (!result.is_valid()) [[unlikely]] {
        return found;
      }
    }

    for (std::size_t i = 0u; i < n_sol; ++i) {
      // Resolve the masks that belong to the surface
      for (const auto &mask :
           detray::ranges::subrange(mask_group, mask_range)) {
        intersection_t is{};

        // Build the resulting intersection(s) from the intersection point
        if constexpr (n_sol > 1) {
          resolve_mask(is, traj, result[i], sf_desc, mask, ctf, cfg,
                       external_mask_tolerance);
        } else {
          resolve_mask(is, traj, result, sf_desc, mask, ctf, cfg,
                       external_mask_tolerance);
        }

        if (is.is_probably_inside()) {
          found[i] = is;
          break;
        }
      }
    }

    return found;
  }
};

template <typename intersection_t, typename... allocator_t>
DETRAY_HOST_DEVICE void insert_sorted(
    const intersection_t &sfi,
    std::vector<intersection_t, allocator_t...> &intersections) {
  auto itr_pos =
      detray::upper_bound(intersections.cbegin(), intersections.cend(), sfi);

  intersections.insert(itr_pos, sfi);
}

/// Specialization for the navigation state cache
template <typename nav_state_t>
DETRAY_HOST_DEVICE void insert_sorted(
    const typename nav_state_t::value_type &sfi, nav_state_t &intersections) {
  auto itr_pos{intersections.cbegin()};

  // For just two candidates int the cache, the navigation state keeps
  // the first as the previously visited candidate -> no sorting needed
  if constexpr (nav_state_t::capacity() > 2u) {
    itr_pos =
        detray::upper_bound(intersections.cbegin(), intersections.cend(), sfi);
  }

  intersections.insert(itr_pos, sfi);
}

/// Add all valid intersections between the trajectory @param traj and the
/// surface @param sf_desc to the container @param is_container .
///
/// Resolves the mask type of the surface in the mask store and runs
/// @c intersection_initialize_impl on the mask group that was found. The
/// resulting candidates are sorted into the container here, so that the
/// container code is instantiated once, rather than once per mask type.
///
/// @tparam intersector_t how to intersect the surface (e.g. ray or helix)
///
/// @param mask_store the store that holds the masks of all surfaces
/// @param mask_link the link to the mask(s) of the surface
/// @param is_container the intersection container to be filled
/// @param traj the input trajectory
/// @param sf_desc the descriptor of the surface to be intersected
/// @param contextual_transforms the store that holds the surface transforms
/// @param ctx the geometry context
/// @param cfg the intersection configuration
/// @param external_mask_tolerance additional mask tol. given by the caller
template <template <typename, typename, bool> class intersector_t,
          typename mask_container_t, typename mask_link_t,
          typename is_container_t, typename traj_t, typename surface_t,
          typename transform_container_t, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE inline void intersection_initialize(
    const mask_container_t &mask_store, const mask_link_t &mask_link,
    is_container_t &is_container, const traj_t &traj, const surface_t &sf_desc,
    const transform_container_t &contextual_transforms,
    const typename transform_container_t::context_type &ctx,
    const intersection::config &cfg,
    const scalar_t external_mask_tolerance = 0.f) {
  using intersection_t = typename is_container_t::value_type;

  const auto found = mask_store.template visit<
      intersection_initialize_impl<intersector_t, intersection_t>>(
      mask_link, traj, sf_desc, contextual_transforms, ctx, cfg,
      external_mask_tolerance);

  for (const auto &is : found) {
    if (is.is_probably_inside()) {
      insert_sorted(is, is_container);
    }
  }
}

/// A functor to update the closest intersection between the trajectory and
/// surface
template <template <typename, typename, bool> class intersector_t>
struct intersection_update {
  /// Operator function to update the intersection
  ///
  /// @tparam mask_group_t is the input mask group type found by variadic
  /// unrolling
  /// @tparam traj_t is the input trajectory type (e.g. ray or helix)
  /// @tparam surface_t is the input surface type
  /// @tparam transform_container_t is the input transform store type
  ///
  /// @param mask_group is the input mask group
  /// @param mask_range is the range of masks in the group that belong to the
  ///                   surface
  /// @param traj is the input trajectory
  /// @param surface is the input surface
  /// @param contextual_transforms is the input transform container
  /// @param mask_tolerance is the tolerance for mask size
  /// @param overstep_tol negative cutoff for the path
  ///
  /// @return the intersection
  template <typename mask_group_t, typename mask_range_t, typename traj_t,
            typename intersection_t, typename transform_container_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline bool operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      const traj_t &traj, intersection_t &sfi,
      const transform_container_t &contextual_transforms,
      const typename transform_container_t::context_type &ctx,
      const intersection::config &cfg,
      const scalar_t external_mask_tolerance = 0.f) const {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;
    using algebra_t = typename mask_t::algebra_type;

    // Find the point of intersection with the underlying geometry
    const auto &ctf = contextual_transforms.at(sfi.surface().transform(), ctx);

    constexpr intersector_t<shape_t, algebra_t, intersection_t::contains_pos()>
        intersector{};
    constexpr std::uint8_t n_sol{decltype(intersector)::n_solutions};

    typename decltype(intersector)::result_type result{};

    if constexpr (concepts::cylindrical<mask_t>) {
      dindex mask_idx{detail::invalid_value<dindex>()};
      if constexpr (concepts::interval<mask_range_t>) {
        mask_idx = mask_range.lower();
      } else {
        mask_idx = mask_range;
      }
      assert(mask_idx < mask_group.size());

      result = intersector.point_of_intersection(
          traj, ctf, mask_group[mask_idx], cfg.overstep_tolerance);
    } else {
      result =
          intersector.point_of_intersection(traj, ctf, cfg.overstep_tolerance);
    }

    // Check if any valid solutions were found
    if constexpr (n_sol > 1) {
      bool found_any{false};
      for (const auto &ip : result) {
        if (ip.is_valid()) {
          found_any = true;
        }
      }
      if (!found_any) [[unlikely]] {
        return false;
      }
    } else {
      if (!result.is_valid()) [[unlikely]] {
        return false;
      }
    }

    // Run over the masks that belong to the surface
    for (const auto &mask : detray::ranges::subrange(mask_group, mask_range)) {
      // Build the resulting intersecion(s) from the intersection point
      if constexpr (n_sol > 1) {
        resolve_mask(sfi, traj, result[0], sfi.surface(), mask, ctf, cfg,
                     external_mask_tolerance);
      } else {
        resolve_mask(sfi, traj, result, sfi.surface(), mask, ctf, cfg,
                     external_mask_tolerance);
      }

      if (sfi.is_probably_inside()) {
        return true;
      }
    }

    return false;
  }
};

}  // namespace detray::detail
