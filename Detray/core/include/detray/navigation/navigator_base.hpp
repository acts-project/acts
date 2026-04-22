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
  /// @param resolve_overstepping whether to check for overstepping
  template <typename track_t, typename nav_state_t, typename context_t>
  DETRAY_HOST_DEVICE constexpr void init(
      const track_t &track, nav_state_t &navigation,
      const navigation::config &cfg, const context_t &ctx,
      const bool resolve_overstepping = false) const {
    // Run local navigation in the current volume
    navigation::local_navigation(track, navigation, cfg, ctx,
                                 resolve_overstepping);

    DETRAY_VERBOSE_HOST("Status: " << navigation.status() << " (next sf.: "
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
  }

  /// @brief Complete update of the navigation flow.
  ///
  /// Restores 'full trust' state to the candidates cache and checks whether
  /// the track stepped onto a portal and a volume switch is due. If so, or
  /// when the previous update according to the given trust level
  /// failed to restore trust, it performs a complete reinitialization of the
  /// navigation.
  ///
  /// @tparam track_t type of track, needs to provide pos() and dir() methods
  ///
  /// @param track access to the track parameters
  /// @param state the current navigation state
  /// @param cfg the navigation configuration
  /// @param ctx the geometry context
  ///
  /// @returns a heartbeat to indicate if the navigation is still alive
  template <typename track_t, typename nav_state_t, typename context_t>
  DETRAY_HOST_DEVICE DETRAY_INLINE constexpr bool update(
      const track_t &track, nav_state_t &navigation,
      const navigation::config &cfg, const context_t &ctx) const {
    assert(navigation.is_alive());
    assert(!track.is_invalid());

    // Update was completely successful (most likely case)
    if (navigation.trust_level() == navigation::trust_level::e_full) {
      DETRAY_VERBOSE_HOST_DEVICE(
          "-> Full trust, nothing left to do: dist to next %f mm",
          navigation());
      return false;
    }

    // Candidates are re-evaluated based on the current trust level.
    // Should result in 'full trust'
    constexpr const navigator_impl_t navigation_impl{};
    bool is_init = navigation_impl.update_impl(track, navigation, cfg, ctx);

    // Otherwise: if we encountered a portal, perform volume switch
    if (navigation.is_on_portal()) {
      navigation::volume_switch(track, navigation, cfg, ctx);
      is_init = true;

      // Reached end of detector
      if (!navigation.is_alive()) {
        return false;
      }
    }
    // If no trust could be restored during the update, try to rescue the
    // navigation stream by re-initializing with loose tolerances
    if (navigation.trust_level() != navigation::trust_level::e_full ||
        navigation.cache_exhausted()) {
      is_init = true;
      navigation::init_loose_cfg(track, navigation, cfg, ctx);

      navigation.run_inspector(cfg, track.pos(), track.dir(), "Re-init: ");
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

    return is_init;
  }
};

}  // namespace detray
