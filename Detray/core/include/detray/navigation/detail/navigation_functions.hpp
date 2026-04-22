// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/logging.hpp"

namespace detray::navigation {

/// A functor that fills the navigation candidates cache by intersecting
/// a surface in the track position neighbourhood with the tangential to
/// the track direction
struct candidate_search {
  /// Visit the volume acceleration data structure(s): Returns a range of
  /// surfaces each which is passed to this functor
  ///
  /// @param sf_descr descriptor of one of the surfaces in the track
  ///                 neighbourhood
  /// @param det access to the detector
  /// @param ctx the geometry context
  /// @param track the track parameters
  /// @param nav_state state of navigation stream of the track
  /// @param mask_tol min. and max. mask tolerance
  /// @param mask_tol_scalor scale factor with track distance in range of
  ///                        @c mask_tol
  /// @param overstep_tol how far behind the track pos to look for
  /// candidates
  template <typename track_t, typename detector_t, typename navigation_state_t>
  DETRAY_HOST_DEVICE constexpr void operator()(
      const typename detector_t::surface_type &sf_descr, const detector_t &det,
      const typename detector_t::geometry_context &ctx, const track_t &track,
      navigation_state_t &nav_state,
      const intersection::config &inter_cfg) const {
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;

    const auto sf = detray::geometry::surface{det, sf_descr};

    DETRAY_DEBUG_HOST("--> Testing surface:\n" << sf);

    // Tangential to the track direction
    detray::detail::ray<algebra_t> tangential{
        track.pos(),
        static_cast<scalar_t>(nav_state.direction()) * track.dir()};

    // Perform intersection and add result to the navigation cache via
    // @c nav_state.insert()
    sf.template visit_mask<
        detray::detail::intersection_initialize<ray_intersector>>(
        nav_state, tangential, sf_descr, det.transform_store(), ctx, inter_cfg,
        nav_state.external_tol());
  }

  /// Test the volume links
  template <typename track_t, typename detector_t, typename navigation_state_t>
  DETRAY_HOST_DEVICE void operator()(
      const dindex & /*vol_idx*/, const detector_t & /*det*/,
      const typename detector_t::geometry_context & /*ctx*/,
      const track_t & /*track*/, navigation_state_t & /*nav_state*/,
      const intersection::config & /*inter_cfg*/) const {
    // Do not search for daughter volumes
  }
};

/// @brief Helper method that updates the intersection of a single candidate
/// and checks reachability
///
/// @tparam candidate_t type of navigation candidate (intersection result)
/// @tparam track_t type of track, needs to provide pos() and dir() methods
/// @tparam detector_t type of the detector
///
/// @param nav_dir the navigation direction (forward/backward)
/// @param candidate the candidate intersection to be updated
/// @param track access to the track parameters
/// @param det access to the detector (geometry)
/// @param cfg the navigation configuration
/// @param external_mask_tolerance additional mask tol. given by the caller
/// @param ctx the geometry context
///
/// @returns @c true if the track can reach this candidate.
template <typename candidate_t, typename track_t, typename detector_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr bool update_candidate(
    const navigation::direction nav_dir, candidate_t &candidate,
    const track_t &track, const detector_t &det,
    const intersection::config &cfg,
    const typename detector_t::scalar_type external_mask_tolerance,
    const typename detector_t::geometry_context &ctx) {
  DETRAY_VERBOSE_HOST_DEVICE("-> Updating target/candidate surface...");

  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;

  // Invalid intersection result cannot be updated
  if (candidate.surface().identifier().is_invalid()) [[unlikely]] {
    return false;
  }

  const auto sf = detray::geometry::surface{det, candidate.surface()};

  // Tangential to the track direction
  auto tangential{detray::detail::ray<algebra_t>(
      track.pos(), static_cast<scalar_t>(nav_dir) * track.dir())};

  // Perform intersection and check whether this candidate is reachable by
  // the track
  return sf.template visit_mask<
      detray::detail::intersection_update<ray_intersector>>(
      std::move(tangential), candidate, det.transform_store(), ctx, cfg,
      external_mask_tolerance);
}

/// @returns @c true if the candidate lies on a surface
template <typename candidate_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr bool has_reached_candidate(
    const candidate_t &candidate, const navigation::config &cfg) {
  return (math::fabs(candidate.path()) < cfg.intersection.path_tolerance);
}

/// @brief Helper method that re-establishes the navigation status after an
/// update.
///
/// It checks weather the track has reached a surface or is still moving
/// towards the next surface candidate. If no new next candidate can be
//  found, it flags 'no trust' in order to trigger a volume initialization.
///
/// @param navigation the current navigation state
/// @param cfg the navigation configuration
template <typename navigation_state_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr void update_status(
    navigation_state_t &navigation, const navigation::config &cfg) {
  DETRAY_VERBOSE_HOST_DEVICE("-> Updating navigation status...");

  // Check whether the track reached the current candidate. Might be a
  // portal, in which case the navigation needs to be re-initialized
  if (!navigation.cache_exhausted() &&
      navigation::has_reached_candidate(navigation.target(), cfg)) {
    navigation.status((navigation.target().surface().is_portal())
                          ? navigation::status::e_on_portal
                          : navigation::status::e_on_object);
    // Set the next object that we want to reach (this function is only
    // called once the state has been freshly updated).
    // Might lead to exhausted cache!
    navigation.advance();
  } else {
    // Otherwise we are moving towards the target
    navigation.status(navigation::status::e_towards_object);
  }
  // Exhaustion happens when after an update no next candidate in the
  // cache is reachable anymore -> triggers init of [new] volume
  // Note: In backwards navigation or with strongly bent tracks, the cache may
  // not be exhausted when trying to exit the volume (the ray is seeing
  // the opposite side of the volume)
  navigation.trust_level(navigation.cache_exhausted() ||
                                 navigation.is_on_portal()
                             ? navigation::trust_level::e_no_trust
                             : navigation::trust_level::e_full);
}

/// @brief Helper method to initialize a navigation state in a given volume.
///
/// Calls the volumes acceleration structure, then tests the surfaces for
/// intersection and keeps the clostest one(s) ("local navigation" in the
/// volume). The closest candidate is set as 'next candidate' or 'target'.
///
/// @tparam track_t type of track, needs to provide pos() and dir() methods
/// @tparam navigation_state_t the state type of the navigation stream
/// @tparam context_t the type of geometry context
///
/// @param track access to the track parameters
/// @param navigation the current navigation state
/// @param cfg the navigation configuration
/// @param ctx the geometry context
template <typename track_t, typename navigation_state_t, typename context_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr void local_navigation(
    const track_t &track, navigation_state_t &navigation,
    const navigation::config &cfg, const context_t &ctx,
    const bool resolve_overstepping = true) {
  DETRAY_VERBOSE_HOST_DEVICE("-> (Re-)initialize detector volume (idx: %d)",
                             navigation.volume());

  // Do not resurrect a failed/finished navigation state
  assert(navigation.is_alive() ||
         navigation.status() == navigation::status::e_unknown);
  assert(!track.is_invalid());

  const auto &det = navigation.detector();
  const auto volume = detray::tracking_volume{det, navigation.volume()};

  // Clean up state
  navigation.clear_cache();

  // Overstepping resolution needed? (if not, make sure full tolerance band is
  // observed around surface)
  intersection::config intr_cfg{cfg.intersection};
  intr_cfg.overstep_tolerance = resolve_overstepping
                                    ? cfg.intersection.overstep_tolerance
                                    : -cfg.intersection.path_tolerance;

  // Search for neighboring surfaces and fill candidates into cache
  using volume_t = typename std::remove_cvref_t<decltype(det)>::volume_type;
  volume.template visit_neighborhood<volume_t::object_id::e_all,
                                     candidate_search>(
      track, cfg.search_window, ctx, det, ctx, track, navigation, intr_cfg);

  // Determine overall state of the navigation after updating the cache
  navigation::update_status(navigation, cfg);

  // If not successful, the propagation setup might be broken
  if (navigation.trust_level() != navigation::trust_level::e_full)
      [[unlikely]] {
    // Do not exit if backward navigation starts on the outmost portal
    if (navigation.is_on_portal()) {
      DETRAY_DEBUG_HOST_DEVICE(
          "-> Adjust trust lvl for possible \"end-of-world\"...");
      navigation.trust_level(
          detray::detail::is_invalid_value(navigation.current().volume_link())
              ? navigation::trust_level::e_full
              : navigation::trust_level::e_no_trust);
    } else if (!navigation.is_on_portal()) {
      DETRAY_VERBOSE_HOST_DEVICE("-> Unable to initialize state!");
    }
  }

  navigation.run_inspector(cfg, track.pos(), track.dir(), "Init complete: ");
}

/// @brief Perform a detector volume switch.
///
/// Once a portal is reached, this function will update the navigation
/// stream to continue in the new volume. If it is a valid detector volume,
/// the navigation is re-initialized by performing local navigation in the
/// volume. If the end of the detector geometry was reached, the navigation
/// exits.
///
/// @tparam track_t type of track, needs to provide pos() and dir() methods
/// @tparam navigation_state_t the state type of the navigation stream
/// @tparam context_t the type of geometry context
///
/// @param track access to the track parameters
/// @param navigation the current navigation state
/// @param cfg the navigation configuration
/// @param ctx the geometry context
template <typename track_t, typename navigation_state_t, typename context_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr void volume_switch(
    const track_t &track, navigation_state_t &navigation,
    const navigation::config &cfg, const context_t &ctx) {
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

  // Initialize new volume. Still on portal: No need to observe overstepping
  local_navigation(track, navigation, cfg, ctx);

  // Fresh initialization, reset trust even though we are on [inner] portal
  navigation.trust_level(navigation::trust_level::e_full);

  DETRAY_VERBOSE_HOST_DEVICE("-> Switched to volume %d", navigation.volume());
}

/// @brief Initialize the volume with loose configuration.
///
/// If trust cannot be established and/or no surfaces can be found in the
/// current volume anymore, try to save the navigation stream by looking
/// for candidates further behind the track position.
///
/// @tparam track_t type of track, needs to provide pos() and dir() methods
/// @tparam navigation_state_t the state type of the navigation stream
/// @tparam context_t the type of geometry context
///
/// @param track access to the track parameters
/// @param navigation the current navigation state
/// @param loose_cfg the navigation configuration (copy on function stack)
/// @param ctx the geometry context
template <typename track_t, typename navigation_state_t, typename context_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr void init_loose_cfg(
    const track_t &track, navigation_state_t &navigation,
    navigation::config loose_cfg, const context_t &ctx) {
  if (navigation.trust_level() != navigation::trust_level::e_full) {
    DETRAY_VERBOSE_HOST_DEVICE("Full trust could not be restored!");
  } else if (navigation.cache_exhausted()) {
    DETRAY_VERBOSE_HOST_DEVICE("Cache exhausted!");
  }
  DETRAY_VERBOSE_HOST_DEVICE("RESCURE MODE: Run init with large tolerances");

  // Use the max mask tolerance in case a track leaves the volume
  // when a sf is 'sticking' out of the portals due to the tol
  const auto new_overstep_tol{
      math::min(100.f * loose_cfg.intersection.overstep_tolerance,
                -10.f * loose_cfg.intersection.max_mask_tolerance)};
  loose_cfg.intersection.overstep_tolerance = new_overstep_tol;

  constexpr bool resolve_overstepping{true};
  local_navigation(track, navigation, loose_cfg, ctx, resolve_overstepping);

  // Unrecoverable
  if (navigation.trust_level() != navigation::trust_level::e_full ||
      navigation.cache_exhausted()) [[unlikely]] {
    navigation.abort("No reachable surfaces");
  }
}

}  // namespace detray::navigation
