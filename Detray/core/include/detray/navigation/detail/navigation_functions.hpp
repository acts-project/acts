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
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/logging.hpp"

namespace detray::navigation {

/// @brief Check whether a candidate is the portal of a daughter volume.
///
/// A volume finds its daughter volumes through their outer portals. These
/// portals belong to the daughter volume and their volume links point away
/// from the daughter volume. Seen from the current (mother) volume, the
/// portal leads into the daughter volume instead.
///
/// @param current_volume index of the volume the navigation is in
/// @param candidate the candidate intersection
///
/// @returns @c true if the candidate surface belongs to a different volume
template <typename candidate_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr bool is_daughter_entry(
    const dindex current_volume, const candidate_t &candidate) {
  return static_cast<dindex>(candidate.surface().volume()) != current_volume;
}

/// @brief Check whether a surface is a portal of a daughter volume.
///
/// The portals of daughter volumes use the mask tolerances (like sensitive
/// surfaces): The tolerance scales with the distance to the portal, so that
/// a daughter volume is not lost while the track is far away, but the portal
/// check is precise when the track reaches the portal. The portals of
/// top-level volumes are checked without tolerance.
///
/// @param det the detector
/// @param sf_descr the surface descriptor
///
/// @returns @c true if the surface is a portal of a volume with a mother
template <concepts::detector detector_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr bool is_daughter_volume_portal(
    const detector_t &det, const typename detector_t::surface_type &sf_descr) {
  return sf_descr.is_portal() && det.volume(sf_descr.volume()).has_mother();
}

/// @param current_volume index of the volume the navigation is in
/// @param candidate a portal candidate
///
/// @returns the index of the volume that the track enters through the portal
/// (invalid, if the track leaves the detector)
template <typename candidate_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr dindex next_volume_index(
    const dindex current_volume, const candidate_t &candidate) {
  if (is_daughter_entry(current_volume, candidate)) {
    return static_cast<dindex>(candidate.surface().volume());
  }
  const auto vol_link{candidate.volume_link()};
  return detray::detail::is_invalid_value(vol_link)
             ? dindex_invalid
             : static_cast<dindex>(vol_link);
}

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
  /// @param tangential the tangential to the track direction
  /// @param nav_state state of navigation stream of the track
  /// @param mask_tol min. and max. mask tolerance
  /// @param mask_tol_scalor scale factor with track distance in range of
  ///                        @c mask_tol
  /// @param overstep_tol how far behind the track pos to look for
  /// candidates
  template <typename traj_t, concepts::detector detector_t,
            typename navigation_state_t>
  DETRAY_HOST_DEVICE constexpr void operator()(
      const typename detector_t::surface_type &sf_descr, const detector_t &det,
      const typename detector_t::geometry_context &ctx,
      const traj_t &tangential, navigation_state_t &nav_state,
      const intersection::config &inter_cfg) const {
    const auto sf = detray::geometry::surface{det, sf_descr};

    DETRAY_DEBUG_HOST("--> Testing surface:\n" << sf);

    // Perform intersection and add result to the navigation cache via
    // @c nav_state.insert()
    sf.template visit_mask<
        detray::detail::intersection_initialize<ray_intersector>>(
        nav_state, tangential, sf_descr, det.transform_store(), ctx, inter_cfg,
        nav_state.external_tol(),
        navigation::is_daughter_volume_portal(det, sf_descr));
  }

  /// Visit a daughter volume: Intersect the outer portals of the daughter
  /// volume. The resulting candidates lead into the daughter volume
  /// (@see @c is_daughter_entry )
  ///
  /// @param daughter_idx index of a daughter volume in the track position
  ///                     neighbourhood
  /// @param det access to the detector
  /// @param ctx the geometry context
  /// @param tangential the tangential to the track direction
  /// @param nav_state state of navigation stream of the track
  /// @param inter_cfg the intersection configuration
  template <typename traj_t, concepts::detector detector_t,
            typename navigation_state_t>
  DETRAY_HOST_DEVICE void operator()(
      const dindex &daughter_idx, const detector_t &det,
      const typename detector_t::geometry_context &ctx,
      const traj_t &tangential, navigation_state_t &nav_state,
      const intersection::config &inter_cfg) const {
    DETRAY_DEBUG_HOST("--> Testing daughter volume: " << daughter_idx);

    const auto daughter = detray::tracking_volume{det, daughter_idx};
    assert(daughter.mother() == nav_state.volume());

    for (const auto &pt_descr : daughter.portals()) {
      (*this)(pt_descr, det, ctx, tangential, nav_state, inter_cfg);
    }
  }
};

/// @brief Check whether the track moves into a daughter volume.
///
/// The navigation can miss the entry portal of a daughter volume, e.g. if the
/// track only clips a corner of the daughter volume. It then reaches the exit
/// portal of the daughter volume, which looks like an entry portal from the
/// mother volume. Check the position slightly ahead of the track.
///
/// @param det the detector
/// @param daughter_idx index of the daughter volume
/// @param track the track (position and direction)
/// @param nav_dir the navigation direction
/// @param ctx the geometry context
///
/// @returns @c false only if the position ahead of the track is outside of
/// the daughter volume
template <concepts::detector detector_t, typename track_t>
DETRAY_HOST_DEVICE constexpr bool is_entering_volume(
    const detector_t &det, const dindex daughter_idx, const track_t &track,
    const navigation::direction nav_dir,
    const typename detector_t::geometry_context &ctx) {
  using scalar_t = typename detector_t::scalar_type;

  // Distance ahead of the track that is tested
  constexpr scalar_t test_dist{10.f * unit<scalar_t>::um};
  constexpr scalar_t tol{1.f * unit<scalar_t>::um};

  const auto ahead =
      track.pos() + (static_cast<scalar_t>(nav_dir) * test_dist) * track.dir();

  return detray::tracking_volume{det, daughter_idx}.is_inside(
             ahead, ctx, tol) != volume_containment::e_outside;
}

/// @brief Helper method that updates the intersection of a single candidate
/// and checks reachability
///
/// @tparam candidate_t type of navigation candidate (intersection result)
/// @tparam traj_t type of trajectory that is intersected with the surface
/// @tparam detector_t type of the detector
///
/// @param candidate the candidate intersection to be updated
/// @param tangential the tangential to the track direction
/// @param det access to the detector (geometry)
/// @param cfg the navigation configuration
/// @param external_mask_tolerance additional mask tol. given by the caller
/// @param ctx the geometry context
///
/// @returns @c true if the track can reach this candidate.
template <typename candidate_t, typename traj_t, concepts::detector detector_t>
DETRAY_HOST_DEVICE DETRAY_INLINE constexpr bool update_candidate(
    candidate_t &candidate, const traj_t &tangential, const detector_t &det,
    const intersection::config &cfg,
    const typename detector_t::scalar_type external_mask_tolerance,
    const typename detector_t::geometry_context &ctx) {
  DETRAY_VERBOSE_HOST_DEVICE("-> Updating target/candidate surface...");

  // Invalid intersection result cannot be updated
  if (candidate.surface().identifier().is_invalid()) [[unlikely]] {
    return false;
  }

  const auto sf = detray::geometry::surface{det, candidate.surface()};

  // Perform intersection and check whether this candidate is reachable by
  // the track
  return sf.template visit_mask<
      detray::detail::intersection_update<ray_intersector>>(
      tangential, candidate, det.transform_store(), ctx, cfg,
      external_mask_tolerance,
      navigation::is_daughter_volume_portal(det, candidate.surface()));
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
DETRAY_HOST_DEVICE
// Not inlining this function gives NVIDIA GPUs more opportunities to
// reconverge threads when this (very large) function is called.
#if defined(__CUDACC__) || defined(__HIP__)
    DETRAY_NO_INLINE
#else
    DETRAY_INLINE
#endif
    constexpr void local_navigation(const track_t &track,
                                    navigation_state_t &navigation,
                                    const navigation::config &cfg,
                                    const context_t &ctx,
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

  using detector_t = typename navigation_state_t::detector_type;
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;

  detray::detail::ray<algebra_t> tangential{
      track.pos(), static_cast<scalar_t>(navigation.direction()) * track.dir()};

  // Search for neighboring surfaces and fill candidates into cache
  using volume_t = typename detector_t::volume_type;
  volume.template visit_neighborhood<volume_t::object_id::e_all,
                                     candidate_search>(
      track, cfg.search_window, ctx, det, ctx, tangential, navigation,
      intr_cfg);

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
          detray::detail::is_invalid_value(navigation::next_volume_index(
              navigation.volume(), navigation.current()))
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
  // The next volume: Either given by the portal link, or the daughter volume
  // that the portal belongs to
  const dindex next_volume{
      navigation::next_volume_index(navigation.volume(), navigation.current())};

  // Navigation reached the end of the detector world
  if (detray::detail::is_invalid_value(next_volume)) [[unlikely]] {
    DETRAY_VERBOSE_HOST_DEVICE("Reached end of detector:");
    navigation.exit();
    return;
  }

  const auto &det = navigation.detector();

  // The volume that the track leaves
  const dindex left_volume{navigation.volume()};

  // Initialize new volume. Still on portal: No need to observe overstepping.
  // This also excludes the portals of a daughter volume, which the track just
  // left, from the candidates (they lie behind the track)
  constexpr bool resolve_overstepping{false};

  // The track reached the portal of a daughter volume, but does not move
  // into it (the entry portal was missed): Stay in the current volume
  if (navigation::is_daughter_entry(left_volume, navigation.current()) &&
      !navigation::is_entering_volume(
          det, next_volume, track, navigation.direction(), ctx)) [[unlikely]] {
    DETRAY_VERBOSE_HOST_DEVICE(
        "-> Track does not enter daughter volume %d: Stay in volume %d",
        next_volume, left_volume);

    local_navigation(track, navigation, cfg, ctx, resolve_overstepping);
    navigation.trust_level(navigation::trust_level::e_full);
    return;
  }

  // Set volume index to the next volume provided by the portal
  navigation.set_volume(next_volume);
  // Check valid volume index
  assert(navigation.volume() < det.volumes().size());

  local_navigation(track, navigation, cfg, ctx, resolve_overstepping);

  // The track can enter a daughter volume at the same position, if the
  // daughter volume touches the boundary of its mother volume: Switch again.
  // Every switch goes one level deeper in the volume hierarchy.
  const auto n_volumes{static_cast<dindex>(det.volumes().size())};

  // The track cannot enter the volume that it leaves, or one of the mother
  // volumes of that volume: It was inside all of them before the switch
  auto is_left_behind = [&det, left_volume, n_volumes](const dindex vol_idx) {
    dindex current{left_volume};
    for (dindex i = 0u; i < n_volumes && current != dindex_invalid; ++i) {
      if (current == vol_idx) {
        return true;
      }
      current = det.volume(current).mother();
    }
    return false;
  };

  for (dindex i = 0u; i < n_volumes && navigation.is_on_portal(); ++i) {
    const auto &candidate = navigation.current();
    const auto daughter_idx{static_cast<dindex>(candidate.surface().volume())};
    if (!navigation::is_daughter_entry(navigation.volume(), candidate) ||
        is_left_behind(daughter_idx) ||
        !navigation::is_entering_volume(det, daughter_idx, track,
                                        navigation.direction(), ctx)) {
      break;
    }

    DETRAY_VERBOSE_HOST_DEVICE("-> Enter daughter volume %d on the portal",
                               daughter_idx);

    navigation.set_volume(daughter_idx);
    local_navigation(track, navigation, cfg, ctx, resolve_overstepping);
  }

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
