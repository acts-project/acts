// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/detail/navigation_functions.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"
#include "detray/navigation/navigator_base.hpp"
#include "detray/utils/logging.hpp"

namespace detray {

namespace navigation {

static constexpr std::size_t default_cache_size{8u};

}  // namespace navigation

/// @brief Navigation class which caches a 'road' through the detector
///
/// This navigator applies a trust level based update of its candidate
/// (intersection) cache, which is kept in the naviagtor's state. The trust
/// level, and with it the appropriate update policy, must be set by an actor,
/// otherwise no update will be performed.
///
/// @tparam detector_t the detector to navigate
/// @tparam k_cache_capacity the capacity of the candidate cache
/// @tparam inspector_t is a validation inspector that can record information
///         about the navigation state at different points of the nav. flow.
/// @tparam intersection_t candidate type
template <typename detector_t,
          std::size_t k_cache_capacity = navigation::default_cache_size,
          typename inspector_t = navigation::void_inspector,
          typename intersection_t = intersection2D<
              typename detector_t::surface_type,
              typename detector_t::algebra_type, !intersection::contains_pos>>
class caching_navigator
    : public navigator_base<caching_navigator<detector_t, k_cache_capacity,
                                              inspector_t, intersection_t>> {
  friend class navigator_base<caching_navigator>;

  using scalar_t = dscalar<typename detector_t::algebra_type>;

 public:
  using detector_type = detector_t;
  using context_type = typename detector_type::geometry_context;
  using intersection_type = intersection_t;
  using inspector_type = inspector_t;

  /// @brief Navigation state that contains a cache of candidates
  ///
  /// Once a volume is reached, the cache is filled by building a 'road' of
  /// surfaces that will likely be encountered by the track in that volume.
  /// The cache keeps a range of reachable candidates that lie between the
  /// next and last index in the cache. The cache is sorted by distance
  /// to the track position.
  class state
      : public navigation::base_state<state, detector_type, k_cache_capacity,
                                      inspector_type, intersection_type> {
    // Allow the navigator classes to access methods that update the state
    friend class caching_navigator;
    friend class navigator_base<caching_navigator>;

    // Allow the filling/updating of candidates
    friend struct detail::intersection_initialize<ray_intersector>;
    friend struct detail::intersection_update<ray_intersector>;

    // Navigation utility functions that need to modify the state
    friend struct navigation::candidate_search;

    template <typename state_t>
    friend constexpr void navigation::update_status(state_t &,
                                                    const navigation::config &);

    template <typename track_t, typename state_t, typename ctx_t>
    friend constexpr void navigation::local_navigation(
        const track_t &, state_t &, const navigation::config &, const ctx_t &,
        const bool);

    template <typename track_t, typename state_t, typename ctx_t>
    friend constexpr void navigation::volume_switch(const track_t &, state_t &,
                                                    const navigation::config &,
                                                    const ctx_t &);

    template <typename track_t, typename state_t, typename ctx_t>
    friend constexpr void navigation::init_loose_cfg(const track_t &, state_t &,
                                                     navigation::config,
                                                     const ctx_t &);
    using base_type =
        navigation::base_state<state, detector_type, k_cache_capacity,
                               inspector_type, intersection_type>;

    // Result of a geometry object evaluation
    using candidate_t = typename base_type::candidate_t;
    using candidate_cache_t = typename base_type::candidate_cache_t;
    using candidate_itr_t = typename base_type::candidate_itr_t;
    using candidate_const_itr_t = typename base_type::candidate_const_itr_t;
    using dist_t = typename base_type::dist_t;

   public:
    using value_type = candidate_t;

    using view_type = detail::get_view_t<inspector_t>;
    using const_view_type = detail::get_view_t<const inspector_t>;

    /// Use common methods of constructing a nvaigation state
    using base_type::base_type;

    /// Update navigation trust level to high trust
    DETRAY_HOST_DEVICE
    constexpr void set_high_trust() {
      DETRAY_VERBOSE_HOST_DEVICE("Flagging re-navigation: \"high trust\"");
      this->trust_level(this->trust_level() < navigation::trust_level::e_high
                            ? this->trust_level()
                            : navigation::trust_level::e_high);
    }

    /// Update navigation trust level to fair trust
    DETRAY_HOST_DEVICE
    constexpr void set_fair_trust() {
      DETRAY_VERBOSE_HOST_DEVICE("Flagging re-navigation: \"fair trust\"");
      this->trust_level(this->trust_level() < navigation::trust_level::e_fair
                            ? this->trust_level()
                            : navigation::trust_level::e_fair);
    }

   private:
    /// Insert a new element @param new_candidate before position @param pos
    DETRAY_HOST_DEVICE
    constexpr void insert(candidate_const_itr_t pos,
                          const intersection_type &new_candidate) {
      // Candidate is too far away to be placed in cache
      if (pos == this->candidates().end()) {
        return;
      }

      assert(detail::is_invalid_value(new_candidate.volume_link()) ||
             new_candidate.volume_link() < this->detector().volumes().size());

      // Insert the first candidate
      if (this->n_candidates() == 0) [[unlikely]] {
        this->candidates()[0] = new_candidate;
        this->last_index(this->last_index() + 1);
        assert(this->next_index() <= this->last_index() + 1);
        assert(static_cast<std::size_t>(this->last_index()) < k_cache_capacity);
        return;
      }

      // Position where to insert the new candidate
      auto idx{static_cast<dist_t>(
          detray::ranges::distance(this->candidates().cbegin(), pos))};
      assert(idx >= 0);

      // Do not add the same surface (intersection) multiple times
      const auto is_overlap_at_pos = [this, &new_candidate](std::size_t index) {
        return (math::fabs(this->candidates()[index].path() -
                           new_candidate.path()) <= 1.f * unit<scalar_t>::um);
      };

      // Do not add the same surface (intersection) multiple times
      const auto is_clash_at_pos = [this, &new_candidate,
                                    &is_overlap_at_pos](std::size_t index) {
        return (this->candidates()[index].surface().identifier() ==
                new_candidate.surface().identifier()) &&
               is_overlap_at_pos(index);
      };

      const auto idxu{static_cast<std::size_t>(idx)};
      if (is_clash_at_pos(idxu) || ((idxu > 0u) && is_clash_at_pos(idxu - 1u)))
          [[unlikely]] {
        return;
      }

      // Shift all following candidates and evict the last element,
      // if the cache is already full
      constexpr auto shift_max{static_cast<dist_t>(k_cache_capacity - 2)};
      const dist_t shift_begin{math::min(this->last_index(), shift_max)};

      // In case of overlaps, prefer sensitives over portals and
      // direct hits over edge hits
      if (is_overlap_at_pos(idxu) &&
          (new_candidate.is_edge() || new_candidate.surface().is_portal())) {
        // Don't shift the overlapping candidate and put the new
        // candidate behind it
        idx++;
        idx = idx > shift_begin ? shift_begin : idx;
      }

      for (dist_t i = shift_begin; i >= idx; --i) {
        const auto j{static_cast<std::size_t>(i)};
        this->candidates()[j + 1u] = this->candidates()[j];
      }

      // Now insert the new candidate and update candidate range
      this->candidates()[static_cast<std::size_t>(idx)] = new_candidate;
      this->last_index(math::min(static_cast<dist_t>(this->last_index() + 1),
                                 static_cast<dist_t>(k_cache_capacity - 1)));

      assert(this->next_index() <= this->last_index() + 1);
      assert(static_cast<std::size_t>(this->last_index()) < k_cache_capacity);
    }

    /// Clear the state
    DETRAY_HOST_DEVICE constexpr void clear_cache() {
      base_type::clear_cache();
      this->next_index(0);
      this->last_index(-1);
    }
  };

 private:
  /// Helper method to update the candidates (surface intersections)
  /// based on an externally provided trust level. Will (re-)initialize the
  /// navigation if there is no trust.
  ///
  /// @tparam track_t type of track, needs to provide pos() and dir() methods
  ///
  /// @param track access to the track parameters
  /// @param state the current navigation state
  /// @param cfg the navigation configuration
  /// @param ctx the geometry context
  template <typename track_t>
  DETRAY_HOST_DEVICE constexpr bool update_impl(const track_t &track,
                                                state &navigation,
                                                const navigation::config &cfg,
                                                const context_type &ctx) const {
    const auto &det = navigation.detector();
    constexpr bool is_init{true};

    assert(navigation.trust_level() != navigation::trust_level::e_full);

    // Update only the current candidate and the corresponding next target
    // - do this only when the navigation state is still coherent
    if (navigation.trust_level() == navigation::trust_level::e_high) {
      DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - high trust");

      // Update next candidate: If not reachable, 'high trust' is broken
      if (!navigation::update_candidate(
              navigation.direction(), navigation.target(), track, det,
              cfg.intersection, navigation.external_tol(), ctx)) {
        DETRAY_VERBOSE_HOST_DEVICE(
            "-> Candidate not reachable! High trust broken:");

        navigation.status(navigation::status::e_unknown);
        // This will run into the fair trust case below.
        navigation.set_fair_trust();
      } else {
        // Update navigation flow on the new candidate information
        navigation::update_status(navigation, cfg);

        navigation.run_inspector(cfg, track.pos(), track.dir(),
                                 "Update complete: high trust: ");

        // The work is done if: the track has not reached a surface yet
        // or trust is gone (portal was reached or the cache is broken).
        if (navigation.status() == navigation::status::e_towards_object ||
            navigation.is_on_portal()) {
          if (navigation.status() == navigation::status::e_towards_object) {
            DETRAY_VERBOSE_HOST_DEVICE("-> Towards surface: %d",
                                       navigation.next_surface().index());
          } else {
            DETRAY_VERBOSE_HOST_DEVICE("-> On portal: idx %d",
                                       navigation.current_surface().index());
          }
          return !is_init;
        }
        // Else (if full trust): Track is on non-portal surface and
        // cache is not exhausted. Ready the next target
        if (navigation.trust_level() == navigation::trust_level::e_full &&
            navigation::update_candidate(
                navigation.direction(), navigation.target(), track, det,
                cfg.intersection, navigation.external_tol(), ctx)) {
          DETRAY_VERBOSE_HOST_DEVICE(
              "-> On non-portal surface (idx %d) and next candidate "
              "in cache is reachable",
              navigation.current_surface().index());
          return !is_init;
        }
        DETRAY_VERBOSE_HOST_DEVICE(
            "-> Next candidate no longer reachable: High trust broken");

        // If next candidate is not reachable, don't 'return', but
        // escalate the trust level.
        // This will run into the fair trust case below or the no trust
        // case if the cache is broken
        navigation.set_fair_trust();
      }
    }
    // Re-evaluate all currently available candidates and sort again
    // - do this when your navigation state is stale, but not invalid
    if (navigation.trust_level() == navigation::trust_level::e_fair &&
        !navigation.cache_exhausted()) {
      DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - fair trust");

      for (auto &candidate : navigation) {
        // Disregard this candidate if it is not reachable
        if (!navigation::update_candidate(navigation.direction(), candidate,
                                          track, det, cfg.intersection,
                                          navigation.external_tol(), ctx)) {
          // Forcefully set dist to numeric max for sorting
          candidate.set_path(std::numeric_limits<scalar_t>::max());
        }
      }
      detray::sequential_sort(navigation.begin(), navigation.end());
      // Take the nearest (sorted) candidate first
      navigation.set_next(navigation.begin());
      // Ignore unreachable elements (needed to determine exhaustion)
      navigation.set_last(find_invalid(navigation.candidates()));
      // Update navigation flow on the new candidate information
      navigation::update_status(navigation, cfg);

      navigation.run_inspector(cfg, track.pos(), track.dir(),
                               "Update complete: fair trust: ");

      // If there are no reachable candidates in the cache after
      // re-evaluation, re-initialize the volume
      if (navigation.cache_exhausted()) {
        navigation.set_no_trust();
      }
    }
    // Re-initialize the volume (actor flagged 'no trust' or previous trust
    // level update failed)
    if (navigation.trust_level() == navigation::trust_level::e_no_trust) {
      DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - no trust");

      constexpr bool resolve_overstepping{true};
      navigation::local_navigation(track, navigation, cfg, ctx,
                                   resolve_overstepping);
      return is_init;
    }

    return !is_init;
  }

  /// Helper to evict all unreachable/invalid candidates from the cache:
  /// Finds the first unreachable candidate (has been invalidated during
  /// update) in a sorted (!) cache.
  ///
  /// @param candidates the cache of candidates to be cleaned
  ///
  /// @returns iterator to the last reachable candidate
  DETRAY_HOST_DEVICE constexpr auto find_invalid(
      const typename state::candidate_cache_t &candidates) const {
    // Depends on previous invalidation of unreachable candidates!
    auto not_reachable = [](const intersection_type &candidate) {
      return candidate.path() == std::numeric_limits<scalar_t>::max();
    };

    return detray::find_if(candidates.begin(), candidates.end(), not_reachable);
  }
};

}  // namespace detray
