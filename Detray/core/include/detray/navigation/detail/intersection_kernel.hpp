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

/// A functor to add all valid intersections between the trajectory and surface
template <template <typename, typename, bool> class intersector_t>
struct intersection_initialize {
  /// Operator function to initialize intersections
  ///
  /// @tparam mask_group_t is the input mask group type found by variadic
  /// unrolling
  /// @tparam is_container_t is the intersection container type
  /// @tparam traj_t is the input trajectory type (e.g. ray or helix)
  /// @tparam surface_t is the input surface type
  /// @tparam transform_container_t is the input transform store type
  ///
  /// @param mask_group is the input mask group
  /// @param is_container is the intersection container to be filled
  /// @param traj is the input trajectory
  /// @param surface is the input surface
  /// @param contextual_transforms is the input transform container
  /// @param mask_tolerance is the tolerance for mask size
  /// @param overstep_tol negative cutoff for the path
  ///
  /// @return the number of valid intersections
  template <typename mask_group_t, typename mask_range_t,
            typename is_container_t, typename traj_t, typename surface_t,
            typename transform_container_t, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline void operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      is_container_t &is_container, const traj_t &traj,
      const surface_t &sf_desc,
      const transform_container_t &contextual_transforms,
      const typename transform_container_t::context_type &ctx,
      const intersection::config &cfg,
      const scalar_t external_mask_tolerance = 0.f) const {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;
    using algebra_t = typename mask_t::algebra_type;
    using intersection_t = typename is_container_t::value_type;

    // Find the point of intersection with the underlying geometry
    const auto &ctf = contextual_transforms.at(sf_desc.transform(), ctx);

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
        return;
      }
    } else {
      if (!result.is_valid()) [[unlikely]] {
        return;
      }
    }

    // Resolve the masks that belong to the surface
    for (const auto &mask : detray::ranges::subrange(mask_group, mask_range)) {
      intersection_t is{};

      // Build the resulting intersecion(s) from the intersection point
      if constexpr (n_sol > 1) {
        std::uint8_t n_found{0u};

        for (std::size_t i = 0u; i < n_sol; ++i) {
          resolve_mask(is, traj, result[i], sf_desc, mask, ctf, cfg,
                       external_mask_tolerance);

          if (is.is_probably_inside()) {
            insert_sorted(is, is_container);
            ++n_found;
          }
          if (n_found == n_sol) {
            return;
          }
        }
      } else {
        resolve_mask(is, traj, result, sf_desc, mask, ctf, cfg,
                     external_mask_tolerance);

        if (is.is_probably_inside()) {
          insert_sorted(is, is_container);
          return;
        }
      }
    }
  }

  template <typename intersection_t>
  DETRAY_HOST_DEVICE void insert_sorted(
      const intersection_t &sfi,
      std::vector<intersection_t> &intersections) const {
    auto itr_pos =
        detray::upper_bound(intersections.cbegin(), intersections.cend(), sfi);

    intersections.insert(itr_pos, sfi);
  }

  /// Specialization for the navigation state cache
  template <typename nav_state_t>
  DETRAY_HOST_DEVICE void insert_sorted(
      const typename nav_state_t::value_type &sfi,
      nav_state_t &intersections) const {
    auto itr_pos{intersections.cbegin()};

    // For just two candidates int the cache, the navigation state keeps
    // the first as the previously visited candidate -> no sorting needed
    if constexpr (nav_state_t::capacity() > 2u) {
      itr_pos = detray::upper_bound(intersections.cbegin(),
                                    intersections.cend(), sfi);
    }

    intersections.insert(itr_pos, sfi);
  }
};

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
