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
#include "detray/geometry/surface.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/detail/navigation_functions.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/navigation/navigation_state.hpp"
#include "detray/utils/logging.hpp"

namespace detray {

/// @brief Navigation base class - interface towards the propagator
///
/// The navigator is initialized around a detector object, but is itself
/// agnostic to the detectors's object/primitive types.
/// Within a detector volume, the navigatior will perform a local navigation
/// based on the geometry acceleration structure(s) that are provided by the
/// volume. Once the local navigation is resolved, it moves to the next volume
/// by a portal.
/// To this end, it requires a link to the [next] navigation volume in every
/// candidate that is computed by intersection from the detector objects:
/// A module surface must link back to its mother volume, while a portal surface
/// links to the next volume in the direction of the track.
///
/// The navigation state is set up by an init() call and then runs an update()
/// after the track state changed.
///
/// @tparam navigator_impl_t implementation of the navigation update
template <typename navigator_impl_t>
class navigator_base {
 public:
  /// @brief Helper method to initialize a volume.
  ///
  /// @tparam track_t type of track, needs to provide pos() and dir() methods
  ///
  /// @param track access to the track parameters
  /// @param state the current navigation state
  /// @param cfg the navigation configuration
  /// @param ctx the geometry context
  ///
  /// @returns whether the navigation is alive after the initialization
  template <typename track_t, typename nav_state_t, typename context_t>
  DETRAY_HOST_DEVICE constexpr bool init(const track_t &track,
                                         nav_state_t &navigation,
                                         const navigation::config &cfg,
                                         const context_t &ctx) const {
    return update(track, navigation, cfg, ctx, navigation::request::e_init);
  }

  /// @brief Complete update of the navigation flow.
  ///
  /// Restores 'full trust' state to the candidates cache and checks whether
  /// the track stepped onto a portal and a volume switch is due. If so, or
  /// when the previous update according to the given trust level
  /// failed to restore trust, it performs a complete reinitialization of the
  /// navigation. Runs all local navigations that @c update_cache requests.
  /// If a local navigation is requested by the caller (@param first), it is
  /// run in place of the cache update, e.g. to initialize the navigation.
  ///
  /// @tparam track_t type of track, needs to provide pos() and dir() methods
  ///
  /// @param track access to the track parameters
  /// @param state the current navigation state
  /// @param cfg the navigation configuration
  /// @param ctx the geometry context
  /// @param first local navigation to run instead of the cache update
  ///
  /// @returns a heartbeat to indicate if the navigation is still alive
  template <typename track_t, typename nav_state_t, typename context_t>
  DETRAY_HOST_DEVICE DETRAY_INLINE constexpr bool update(
      const track_t &track, nav_state_t &navigation,
      const navigation::config &cfg, const context_t &ctx,
      const navigation::request first = navigation::request::e_none) const {
    using enum navigation::request;

    // Restore the trust in the cache, unless a local navigation is due anyway
    navigation::update_result res{first, false};
    if (first == e_none) {
      res = update_cache(track, navigation, cfg, ctx);
    }
    bool is_init = res.is_init;

    // Run the local navigations that are needed to complete the update; if
    // the next request is e_none, then we are done. Otherwise, we determine
    // the next request type from the current one, and the navigation state.
    for (navigation::request req = res.next; req != e_none;
         req = next_request(navigation, req)) {
      perform(req, track, navigation, cfg, ctx);
      is_init = true;

      // Reached end of detector
      if (!navigation.is_alive()) {
        return false;
      }
    }

    return is_init;
  }

  /// @brief Cheap part of the navigation update.
  ///
  /// Restores the trust in the candidates cache as far as possible without a
  /// local navigation and reports the local navigation that is needed to
  /// complete the update, if any. The caller runs it with @c perform and
  /// continues with @c next_request until no request is left.
  ///
  /// @param track access to the track parameters
  /// @param navigation the current navigation state
  /// @param cfg the navigation configuration
  /// @param ctx the geometry context
  ///
  /// @returns the local navigation that has to follow, if any
  template <typename track_t, typename nav_state_t, typename context_t>
  DETRAY_HOST_DEVICE DETRAY_INLINE constexpr navigation::update_result
  update_cache(const track_t &track, nav_state_t &navigation,
               const navigation::config &cfg, const context_t &ctx) const {
    assert(navigation.is_alive());
    assert(!track.is_invalid());

    // Update was completely successful (most likely case)
    if (navigation.trust_level() == navigation::trust_level::e_full) {
      DETRAY_VERBOSE_HOST_DEVICE(
          "-> Full trust, nothing left to do: dist to next %f mm",
          navigation());
      return {};
    }

    // Candidates are re-evaluated based on the current trust level.
    // Should result in 'full trust'
    constexpr const navigator_impl_t navigation_impl{};
    const bool is_init =
        navigation_impl.update_impl(track, navigation, cfg, ctx);

    // Re-initialize the volume (actor flagged 'no trust' or previous trust
    // level update failed). Note: The track can also be on a portal with
    // no trust, which requires a volume switch, not a re-initialization
    if (is_init) {
      DETRAY_VERBOSE_HOST_DEVICE("Called 'update()' - no trust");
      return {navigation::request::e_re_init, false};
    }

    // Otherwise: check whether a volume switch or a rescue is needed
    return {next_request(navigation, navigation::request::e_none), false};
  }

  /// @brief Run the local navigation @param req
  ///
  /// @param req the requested local navigation
  /// @param track access to the track parameters
  /// @param navigation the current navigation state
  /// @param cfg the navigation configuration
  /// @param ctx the geometry context
  template <typename track_t, typename nav_state_t, typename context_t>
  DETRAY_HOST_DEVICE DETRAY_INLINE constexpr void perform(
      const navigation::request req, const track_t &track,
      nav_state_t &navigation, const navigation::config &cfg,
      const context_t &ctx) const {
    using enum navigation::request;
    assert(req != e_none);

    // Rescue mode: run the local navigation with loose tolerances. Use the
    // max mask tolerance in case a track leaves the volume when a sf is
    // 'sticking' out of the portals due to the tol
    navigation::config loose_cfg{};
    const navigation::config *nav_cfg{&cfg};
    if (req == e_loose_re_init) {
      if (navigation.trust_level() != navigation::trust_level::e_full) {
        DETRAY_VERBOSE_HOST_DEVICE("Full trust could not be restored!");
      } else if (navigation.cache_exhausted()) {
        DETRAY_VERBOSE_HOST_DEVICE("Cache exhausted!");
      }
      DETRAY_VERBOSE_HOST_DEVICE(
          "RESCURE MODE: Run init with large tolerances");

      loose_cfg = cfg;
      loose_cfg.intersection.overstep_tolerance =
          math::min(100.f * cfg.intersection.overstep_tolerance,
                    -10.f * cfg.intersection.max_mask_tolerance);
      nav_cfg = &loose_cfg;
    }

    // Volume switch: move to the volume behind the portal
    if (req == e_volume_switch) {
      // Navigation reached the end of the detector world
      if (detray::detail::is_invalid_value(navigation.current().volume_link()))
          [[unlikely]] {
        DETRAY_VERBOSE_HOST_DEVICE("Reached end of detector:");
        navigation.exit();
        return;
      }

      // Set volume index to the next volume provided by the portal
      navigation.set_volume(navigation.current().volume_link());
      // Check valid volume index
      assert(navigation.volume() < navigation.detector().volumes().size());
    }

    // Run the local navigation in the current volume. Overstepping is only
    // resolved once the track is in flight, not at the initial navigation
    const bool resolve_overstepping{req != e_init};
    navigation::local_navigation(track, navigation, *nav_cfg, ctx,
                                 resolve_overstepping);

    if (req == e_volume_switch) {
      // Fresh initialization, reset trust even though we are on [inner] portal
      navigation.trust_level(navigation::trust_level::e_full);

      DETRAY_VERBOSE_HOST_DEVICE("-> Switched to volume %d",
                                 navigation.volume());
    }

    if (req == e_loose_re_init) {
      // Unrecoverable
      if (navigation.trust_level() != navigation::trust_level::e_full ||
          navigation.cache_exhausted()) [[unlikely]] {
        navigation.abort("No reachable surfaces");
      }

      navigation.run_inspector(cfg, track.pos(), track.dir(), "Re-init: ");
    }
  }

  /// @brief Checks whether the track stepped onto a portal and a volume
  /// switch is due, or whether the previous (re-)initialization failed to
  /// restore trust, so that a re-initialization with loose tolerances is
  /// needed.
  ///
  /// @param navigation the current navigation state
  /// @param last the local navigation that was performed last
  ///
  /// @returns the local navigation that has to follow @param last
  template <typename nav_state_t>
  DETRAY_HOST_DEVICE DETRAY_INLINE constexpr navigation::request next_request(
      const nav_state_t &navigation, const navigation::request last) const {
    using enum navigation::request;

    // Nothing follows the end of the detector
    if (!navigation.is_alive()) {
      return e_none;
    }
    // Nothing follows the initial navigation, either
    if (last != e_init) {
      // If we encountered a portal, perform volume switch
      if ((last == e_none || last == e_re_init) && navigation.is_on_portal()) {
        return e_volume_switch;
      }
      // If no trust could be restored during the update, try to rescue the
      // navigation stream by re-initializing with loose tolerances
      if (last != e_loose_re_init &&
          (navigation.trust_level() != navigation::trust_level::e_full ||
           navigation.cache_exhausted())) {
        return e_loose_re_init;
      }
    }

    DETRAY_VERBOSE_HOST("Status: " << navigation.status() << " (vol.:"
                                   << navigation.volume() << ", next sf.: "
                                   << navigation.next_surface().index() << ")");
    if (navigation.is_on_surface()) {
      DETRAY_VERBOSE_HOST("-> Current surface: "
                          << navigation.current_surface().index()
                          << ", has material: " << std::boolalpha
                          << navigation.current_surface().has_material()
                          << std::noboolalpha);
    }
    DETRAY_VERBOSE_HOST_DEVICE("Update complete: dist to next %f mm",
                               navigation());

    return e_none;
  }
};

}  // namespace detray
