// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/detail/navigation_functions.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"
#include "detray/navigation/navigator_base.hpp"
#include "detray/utils/logging.hpp"

namespace detray {

/// @brief The geometry navigation class.
///
/// @tparam detector_t the detector to navigate
/// @tparam inspector_t is a validation inspector that can record information
///         about the navigation state at different points of the nav. flow.
/// @tparam intersection_t candidate type
template <
    typename detector_t, typename inspector_t = navigation::void_inspector,
    typename intersection_t = intersection2D<typename detector_t::surface_type,
                                             typename detector_t::algebra_type,
                                             !intersection::contains_pos>>
class navigator : public navigator_base<
                      navigator<detector_t, inspector_t, intersection_t>> {
  friend class navigator_base<navigator>;

  using scalar_t = dscalar<typename detector_t::algebra_type>;

 public:
  using detector_type = detector_t;
  using context_type = typename detector_type::geometry_context;
  using intersection_type = intersection_t;
  using inspector_type = inspector_t;

  /// @brief Navigation state that contains a the current and the next
  /// candidates
  class state
      : public navigation::base_state<state, detector_type, 2u, inspector_type,
                                      intersection_type> {
    // Allow the navigator classes to access methods that update the state
    friend class navigator;
    friend class navigator_base<navigator>;

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
    using base_type = navigation::base_state<state, detector_type, 2u,
                                             inspector_type, intersection_type>;

    using candidate_const_itr_t = typename base_type::candidate_const_itr_t;

   public:
    using value_type = typename base_type::candidate_t;
    using view_type = detail::get_view_t<inspector_t>;
    using const_view_type = detail::get_view_t<const inspector_t>;

    /// Use common methods of constructing a nvaigation state
    using base_type::base_type;

    /// @returns the navigator always has only one candidate
    DETRAY_HOST_DEVICE
    constexpr auto n_candidates() const -> dindex { return 1u; }

    /// Update navigation trust level to high trust
    /// @{
    DETRAY_HOST_DEVICE
    constexpr void set_high_trust() {
      this->trust_level(this->trust_level() < navigation::trust_level::e_high
                            ? this->trust_level()
                            : navigation::trust_level::e_high);
    }
    DETRAY_HOST_DEVICE
    constexpr void set_fair_trust() { this->set_high_trust(); }
    /// @}

   private:
    /// Insert a new element @param new_cadidate before position @param pos
    DETRAY_HOST_DEVICE
    constexpr void insert(candidate_const_itr_t /*pos*/,
                          const intersection_type &new_cadidate) {
      // Insert the first candidate
      if (math::fabs(new_cadidate) < math::fabs(this->candidates()[1].path)) {
        this->candidates()[1] = new_cadidate;
      }

      assert(this->next_index() == 1u);
      assert(detail::is_invalid_value(new_cadidate.volume_link()) ||
             new_cadidate.volume_link() < this->detector().volumes().size());
    }

    /// If the current target is reached, move it to the first position in
    /// the cache and reset the target slot
    DETRAY_HOST_DEVICE
    constexpr void reset_candidate() {
      // Move the current candidate to the first cache position
      this->candidates()[0] = this->candidates()[1];

      // Flag the old candidate as invalid
      this->candidates()[1].set_path(std::numeric_limits<scalar_t>::max());
    }
  };

 private:
  /// Helper method to update the candidate (surface intersections)
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

    // Update only the current candidate and the corresponding
    // - do this only when the navigation state is still coherent
    if (navigation.trust_level() == navigation::trust_level::e_high) {
      DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - high trust");

      // Update next candidate: If not reachable, 'high trust' is broken
      if (!navigation::update_candidate(
              navigation.direction(), navigation.target(), track, det,
              cfg.intersection, navigation.external_tol(), ctx)) {
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
          return !is_init;
        }
        // Else (if full trust): Track is on non-portal surface and
        // cache is not exhausted. Move the current surface back in the
        // cache and re-initialize the volume to find the next target
        if (navigation.trust_level() == navigation::trust_level::e_full) {
          navigation.reset_candidate();
        }

        // Find the next/different candidate
        navigation.set_no_trust();
      }
    }
    // If [next] target is not reachable or actor flagged 'no trust',
    // re-initialize the volume.
    assert(navigation.trust_level() == navigation::trust_level::e_no_trust);
    DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - no trust");

    constexpr bool resolve_overstepping{true};
    navigation::local_navigation(track, navigation, cfg, ctx,
                                 resolve_overstepping);
    return is_init;
  }
};

}  // namespace detray
