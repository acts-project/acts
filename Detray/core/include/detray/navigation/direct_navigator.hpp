// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/detail/navigation_functions.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"

namespace detray {

template <typename detector_t, typename surface_t = void>
class direct_navigator {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;

 public:
  using detector_type = detector_t;
  using context_type = detector_type::geometry_context;

  // Use an external surface type, if one was provided
  using surface_type =
      std::conditional_t<std::same_as<surface_t, void>,
                         typename detector_t::surface_type, surface_t>;
  using intersection_type =
      intersection2D<surface_type, algebra_t, !intersection::contains_pos>;
  using inspector_type = navigation::void_inspector;

  class state
      : public navigation::base_state<state, detector_type, 2u, inspector_type,
                                      intersection_type> {
    friend class direct_navigator;
    friend struct detail::intersection_update<ray_intersector>;

    template <typename state_t>
    friend constexpr void navigation::update_status(state_t &,
                                                    const navigation::config &);

    using base_type = navigation::base_state<state, detector_type, 2u,
                                             inspector_type, intersection_type>;

    // Internal surface index
    using dist_t = std::int_least16_t;

   public:
    using value_type = intersection_type;
    using sequence_type = vecmem::device_vector<surface_type>;

    using view_type = dvector_view<surface_type>;
    using const_view_type = dvector_view<const surface_type>;

    /// Constructor using the detector and an externally provided sequence
    /// of detector surfaces
    DETRAY_HOST_DEVICE constexpr state(const detector_t &det,
                                       const view_type &sequence)
        : base_type(det), m_sequence(sequence) {
      assert(!m_sequence.empty());
      reset();
    }

    /// Reset the navigation state
    DETRAY_HOST_DEVICE void reset() {
      // Remove old navigation information
      clear_cache();

      // Set the index into the external surface sequence to the beginning
      // or end of the container
      m_next_external =
          is_forward() ? 0 : static_cast<dist_t>(m_sequence.size()) - 1;

      // Update the target with the next external surface
      this->target().set_surface(next_external());

      // The next candidate is always stored in the second cache entry
      this->next_index(1u);
      this->last_index(1u);
    }

    /// @returns the direct navigator always has only one candidate
    DETRAY_HOST_DEVICE
    constexpr auto n_candidates() const -> dindex { return 1u; }

    /// @returns the externally provided mask tolerance - const
    DETRAY_HOST_DEVICE
    constexpr scalar_t external_tol() const {
      return std::numeric_limits<scalar_t>::max();
    }

    /// Set direction in which the navigator should search for candidates
    DETRAY_HOST_DEVICE
    constexpr void set_direction(const navigation::direction dir) {
      base_type::set_direction(dir);
    }

    /// @returns true if the the next external is available in the sequence
    DETRAY_HOST_DEVICE
    constexpr bool has_next_external() const {
      return (is_forward() &&
              m_next_external != static_cast<dist_t>(m_sequence.size())) ||
             (!is_forward() && m_next_external != -1);
    }

    /// Advance the iterator
    DETRAY_HOST_DEVICE
    constexpr surface_type next_external() {
      assert(has_next_external());
      return m_sequence[static_cast<unsigned int>(m_next_external)];
    }

    /// Advance the iterator (navigation status needs to be up to date!)
    DETRAY_HOST_DEVICE
    constexpr void advance() {
      // The target has become the current candidate
      this->candidates()[0] = this->target();

      assert(has_next_external());
      if (is_forward()) {
        m_next_external++;
      } else {
        m_next_external--;
      }

      // Update the target with the next external surface
      // Could make the target invalid -> exit navigation
      if (has_next_external()) {
        // If the next external can be indexed, the cast is safe
        this->target().set_surface(
            m_sequence[static_cast<unsigned int>(m_next_external)]);
      }

      assert(
          !has_next_external() ||
          (has_next_external() && (this->target().surface().is_sensitive() ||
                                   this->target().surface().has_material())));
    }

    /// Clear the state
    DETRAY_HOST_DEVICE constexpr void clear_cache() {
      base_type::clear_cache();
      this->next_index(1);
      this->last_index(1);
    }

    /// @returns flag that indicates whether navigation was successful
    DETRAY_HOST_DEVICE
    constexpr bool finished() const {
      // Normal exit for this navigation?
      bool is_finished = base_type::finished();

      // All external surfaces were visited?
      is_finished &=
          ((is_forward() &&
            m_next_external == static_cast<dist_t>(m_sequence.size())) ||
           (!is_forward() && m_next_external == -1));

      return is_finished;
    }

   private:
    /// @returns 'true' if the navigation direction is 'forward'
    DETRAY_HOST_DEVICE
    constexpr bool is_forward() const {
      return this->direction() == navigation::direction::e_forward;
    }

    /// Target surfaces
    sequence_type m_sequence;

    /// Index of the next (target) surface descriptor in the sequence
    dist_t m_next_external{0};
  };

  /// Initialize the direct navigation flow on the first/last external surface
  template <typename track_t>
  DETRAY_HOST_DEVICE inline void init(const track_t &track, state &navigation,
                                      const navigation::config &cfg,
                                      const context_type &ctx) const {
    DETRAY_VERBOSE_HOST_DEVICE("Called 'init()':");
    assert(navigation.has_next_external());

    // Clean up state
    navigation.clear_cache();

    // Update the next candidate
    update(track, navigation, cfg, ctx);

    DETRAY_VERBOSE_HOST_DEVICE("Init complete!");
  }

  /// Update the navigation status on the current next external and switch
  /// to the following external surface if the current one was reached
  template <typename track_t>
  DETRAY_HOST_DEVICE inline bool update(const track_t &track, state &navigation,
                                        const navigation::config &cfg,
                                        const context_type &ctx = {}) const {
    DETRAY_VERBOSE_HOST_DEVICE("Called 'update()'");
    DETRAY_DEBUG_HOST(" -> Trust level: " << navigation.trust_level());

    // Do not resurrect a failed/finished navigation state
    assert(navigation.is_alive() ||
           navigation.status() == navigation::status::e_unknown);
    assert(!track.is_invalid());

    // If the last external was reached, the navigation is finished
    constexpr bool is_init{true};
    if (!navigation.has_next_external()) {
      DETRAY_VERBOSE_HOST_DEVICE("No next target in surface sequence: Exit");
      navigation.exit();
      return !is_init;
    }

    const detector_type &det = navigation.detector();

    // Update only the current candidate and the corresponding next target
    if (navigation.trust_level() != navigation::trust_level::e_full) {
      // Use infinite tolerance to hit the target surface
      constexpr auto inf{std::numeric_limits<float>::max()};
      const intersection::config intr_cfg{
          inf, inf, inf, cfg.intersection.path_tolerance, -inf};

      // Update the current target. If it cannot be reached, direct
      // navigation is broken
      if (!navigation::update_candidate(
              navigation.direction(), navigation.target(), track, det, intr_cfg,
              navigation.external_tol(), ctx)) {
        navigation.abort("Could not reach current target");
        return !is_init;
      }

      // Update navigation flow on the new candidate information and set
      // the next target if the surface has been reached
      navigation::update_status(navigation, cfg);
      // Set full trust, even if a portal was reached (no volume switch)
      navigation.trust_level(navigation::trust_level::e_full);

      // The work is done if the track has not reached the surface yet
      if (navigation.status() == navigation::status::e_towards_object) {
        DETRAY_VERBOSE_HOST_DEVICE("Update complete (towards object):");
        DETRAY_DEBUG_HOST("\n"
                          << detray::navigation::print_state(navigation)
                          << detray::navigation::print_candidates(
                                 navigation, cfg, track.pos(), track.dir()));

        return !is_init;
      }

      // Otherwise, track is on surface: Update the next target
      if (navigation.has_next_external() &&
          !navigation::update_candidate(
              navigation.direction(), navigation.target(), track, det, intr_cfg,
              navigation.external_tol(), ctx)) {
        navigation.abort("Could not find new target after surface was reached");
        return !is_init;
      }

      // At this point, the track has to be on surface:
      // Set volume index to the volume the current surface is in
      navigation.set_volume(navigation.current().surface().volume());

      // Set full trust again: no volume switch needed
      navigation.trust_level(navigation::trust_level::e_full);

      DETRAY_VERBOSE_HOST_DEVICE("Update complete (on surface):");
      DETRAY_DEBUG_HOST("\n"
                        << detray::navigation::print_state(navigation)
                        << detray::navigation::print_candidates(
                               navigation, cfg, track.pos(), track.dir()));

      // Return true to reset the step size of the RKN algorithm
      return is_init;
    }

    DETRAY_VERBOSE_HOST_DEVICE(" -> Full trust: Nothing left to do");

    return !is_init;
  }
};

}  // namespace detray
