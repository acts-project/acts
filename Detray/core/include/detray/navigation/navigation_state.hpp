// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/navigation.hpp"
#include "detray/geometry/identifier.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/detail/print_state.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <limits>

namespace detray::navigation {

/// @brief A void inpector that does nothing.
///
/// Inspectors can be plugged in to understand the current navigation state.
struct void_inspector {
  struct void_view : public detray::detail::dbase_view {};

  using view_type = void_view;
  using const_view_type = const void_view;

  constexpr void_inspector() = default;

  DETRAY_HOST_DEVICE
  constexpr explicit void_inspector(
      const void_view & /*ignored*/) { /*Do nothing*/ }

  template <typename state_t>
  DETRAY_HOST_DEVICE constexpr void operator()(const state_t & /*ignored*/,
                                               const char * /*ignored*/) const {
    /*Do nothing*/
  }
};

/// @brief A navigation state object used to cache the information of the
/// current navigation stream.
///
/// The state is passed between navigation calls and is accessible to the
/// actors in the propagation, for which it defines the public interface
/// towards the navigation. The navigator class is responsible for updating the
/// elements in the navigation state with every navigation call, establishing
/// 'full trust' after changes to the track state reduced the trust level.
///
/// @tparam detector_t the type of the detector that is being navigated
/// @tparam k_cache_capacity number of object candidates the state can hold
/// @tparam inspector_t is a validation inspector that can record information
///         about the navigation state at different points of the nav. flow.
/// @tparam intersection_t result of an intersection operation
template <typename derived_t, typename detector_t, std::size_t k_cache_capacity,
          typename inspector_t, typename intersection_t>
class base_state {
  /// Need at least two slots for the caching to work
  static_assert(k_cache_capacity >= 2u,
                "Navigation cache needs to have a capacity larger than 1");

  /// Cast to the navigation state type to access its methods
  /// @{
  DETRAY_HOST_DEVICE
  constexpr auto cast_impl() -> derived_t & {
    return static_cast<derived_t &>(*this);
  }
  DETRAY_HOST_DEVICE
  constexpr auto cast_impl() const -> const derived_t & {
    return static_cast<const derived_t &>(*this);
  }
  /// @}

 protected:
  // Linear algebra types
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;

  // Result of a geometry object evaluation
  using candidate_t = intersection_t;
  using candidate_cache_t = darray<candidate_t, k_cache_capacity>;
  using candidate_itr_t = typename candidate_cache_t::iterator;
  using candidate_const_itr_t = typename candidate_cache_t::const_iterator;
  using dist_t = std::int_least8_t;

  // Link type to the next detector volume
  using nav_link_t = typename detector_t::surface_type::navigation_link;

 public:
  using detector_type = detector_t;

  /// Default constructor (needs a detector)
  base_state() = delete;

  /// Constructor using a given detector @param det
  DETRAY_HOST_DEVICE
  constexpr explicit base_state(const detector_t &det) : m_detector(&det) {
    clear_cache();
  }

  /// Construct from detector @param det and inspector view @param view
  template <concepts::device_view view_t>
  DETRAY_HOST_DEVICE constexpr base_state(const detector_t &det, view_t view)
      : m_detector(&det), m_inspector(view) {
    clear_cache();
  }

  /// @returns a reference to the detector
  DETRAY_HOST_DEVICE
  constexpr auto detector() const -> const detector_t & {
    return (*m_detector);
  }

  /// @returns the capacity of the internal candidate storage
  static consteval std::size_t capacity() { return k_cache_capacity; }

  /// @returns current navigation status - const
  DETRAY_HOST_DEVICE
  constexpr auto status() const -> navigation::status { return m_status; }

  /// @returns the navigation heartbeat (is this stream still being updated?)
  DETRAY_HOST_DEVICE
  constexpr bool is_alive() const {
    return cast_impl().status() > navigation::status::e_unknown;
  }

  /// @returns current navigation direction - const
  DETRAY_HOST_DEVICE
  constexpr auto direction() const -> navigation::direction {
    return m_direction;
  }

  /// Set direction in which the navigator should search for candidates
  DETRAY_HOST_DEVICE
  constexpr void set_direction(const navigation::direction dir) {
    if (m_direction != dir) {
      // Force renavigation after direction change
      cast_impl().set_no_trust();
      DETRAY_VERBOSE_HOST_DEVICE("Re-navigate after explicit direction change");

      m_direction = dir;
    }
  }

  /// @returns the externally provided mask tolerance - const
  DETRAY_HOST_DEVICE
  constexpr scalar_t external_tol() const { return m_external_mask_tol; }

  /// Set externally provided mask tolerance according to noise prediction
  DETRAY_HOST_DEVICE
  constexpr void set_external_tol(const scalar_t tol) {
    DETRAY_VERBOSE_HOST("Setting external mask tolerance: " << tol);
    m_external_mask_tol = tol;
  }

  /// @returns navigation trust level - const
  DETRAY_HOST_DEVICE
  constexpr auto trust_level() const -> navigation::trust_level {
    return m_trust_level;
  }

  /// Update navigation trust level to 'no trust'
  /// @{
  DETRAY_HOST_DEVICE
  constexpr void set_no_trust() {
    DETRAY_VERBOSE_HOST_DEVICE("Flagging re-navigation: \"no trust\"");
    m_trust_level = navigation::trust_level::e_no_trust;
  }

  /// @note only two trust levels for basic navigation implementation
  DETRAY_HOST_DEVICE
  constexpr void set_high_trust() { cast_impl().set_no_trust(); }

  /// @note only two trust levels for basic navigation implementation
  DETRAY_HOST_DEVICE
  constexpr void set_fair_trust() { cast_impl().set_no_trust(); }
  /// @}

  /// Helper method to check the track has reached a module surface
  DETRAY_HOST_DEVICE
  constexpr auto is_on_surface() const -> bool {
    return (m_status == navigation::status::e_on_object ||
            m_status == navigation::status::e_on_portal);
  }

  /// Helper method to check the track has reached a sensitive surface
  DETRAY_HOST_DEVICE
  constexpr auto is_on_sensitive() const -> bool {
    return (m_status == navigation::status::e_on_object) &&
           (cast_impl().geometry_identifier().id() == surface_id::e_sensitive);
  }

  /// Helper method to check the track has reached a passive surface
  DETRAY_HOST_DEVICE
  constexpr auto is_on_passive() const -> bool {
    return (m_status == navigation::status::e_on_object) &&
           (cast_impl().geometry_identifier().id() == surface_id::e_passive);
  }

  /// Helper method to check the track has reached a portal surface
  DETRAY_HOST_DEVICE
  constexpr auto is_on_portal() const -> bool {
    return m_status == navigation::status::e_on_portal;
  }

  /// Helper method to check the track has encountered material
  DETRAY_HOST_DEVICE
  constexpr auto encountered_sf_material() const -> bool {
    return (cast_impl().is_on_surface()) &&
           (cast_impl().current().surface().has_material());
  }

  /// @returns flag that indicates whether navigation was successful
  DETRAY_HOST_DEVICE
  constexpr bool finished() const {
    // Normal exit for this navigation?
    return (m_status == navigation::status::e_exit);
  }

  /// @returns current volume index - const
  DETRAY_HOST_DEVICE
  constexpr auto volume() const -> nav_link_t { return m_volume_index; }

  /// Set start/new volume from volume index @param v
  DETRAY_HOST_DEVICE
  constexpr void set_volume(dindex v) {
    assert(detray::detail::is_invalid_value(static_cast<nav_link_t>(v)) ||
           v < cast_impl().detector().volumes().size());
    if (v != m_volume_index) {
      DETRAY_VERBOSE_HOST_DEVICE("Setting new volume %d", v);
      // Make sure the new volume is properly initialized
      cast_impl().set_no_trust();
    }
    m_volume_index = static_cast<nav_link_t>(v);
  }

  /// @returns number of currently cached (reachable) candidates, i.e. the
  /// target and everything behind it - const
  DETRAY_HOST_DEVICE
  constexpr auto n_candidates() const -> dindex {
    assert(m_n_valid >= static_cast<dindex>(m_first));
    return m_n_valid - static_cast<dindex>(m_first);
  }

  /// Set the number of reachable candidates (the current candidate, if any,
  /// stays valid)
  DETRAY_HOST_DEVICE
  constexpr void n_candidates(dindex n) {
    m_n_valid = static_cast<dindex>(m_first) + n;
    assert(m_n_valid <= static_cast<dindex>(usable_capacity()));
  }

  /// @returns true if there are no cached candidates left - const
  DETRAY_HOST_DEVICE
  constexpr bool cache_exhausted() const {
    return cast_impl().n_candidates() == 0u;
  }

  /// @returns current/previous object that was reached - const
  DETRAY_HOST_DEVICE
  constexpr auto current() const -> const candidate_t {
    assert(cast_impl().is_on_surface());
    assert(m_first == 1);
    return this->candidate_at(0);
  }

  /// @returns next object that we want to reach (current target) - const
  DETRAY_HOST_DEVICE
  constexpr auto target() const -> const candidate_t {
    return (m_first == 0) ? this->candidate_at(0) : this->candidate_at(1);
  }

  /// @returns identifier of the detector surface the navigator is on
  /// (invalid when not on surface) - const
  DETRAY_HOST_DEVICE
  constexpr auto geometry_identifier() const -> geometry::identifier {
    return cast_impl().current().surface().identifier();
  }

  /// @returns true if the current candidate lies on the surface edge
  DETRAY_HOST_DEVICE
  constexpr bool is_edge_candidate() const {
    assert(cast_impl().is_on_surface());
    return cast_impl().current().is_edge();
  }

  /// @returns true if the current candidate lies on the surface
  DETRAY_HOST_DEVICE
  constexpr bool is_good_candidate() const {
    assert(cast_impl().is_on_surface());
    return cast_impl().current().is_inside();
  }

  /// @returns true if the current candidate lies on the surface,
  /// including its edge
  DETRAY_HOST_DEVICE
  constexpr bool is_probably_candidate() const {
    assert(cast_impl().is_on_surface());
    return cast_impl().current().is_probably_inside();
  }

  /// Scalar representation of the navigation state,
  /// @returns distance to next target
  DETRAY_HOST_DEVICE
  constexpr scalar_t operator()() const {
    assert(math::isfinite(cast_impl().target().path()));
    return static_cast<scalar_t>(cast_impl().direction()) *
           cast_impl().target().path();
  }

  /// @param surface_t the surface interface type that is required
  /// @returns current detector surface the navigator is on - const
  template <template <typename> class surface_t = tracking_surface>
  DETRAY_HOST_DEVICE constexpr auto current_surface() const {
    assert(cast_impl().is_on_surface());
    return surface_t<detector_t>{*m_detector, current().surface()};
  }

  /// @param volume_t the volume interface type that is required
  /// @returns current detector volume of the navigation stream
  template <template <typename> class volume_t = tracking_volume>
  DETRAY_HOST_DEVICE constexpr auto current_volume() const {
    return volume_t<detector_t>{*m_detector, m_volume_index};
  }

  /// @param surface_t the surface interface type that is required
  /// @returns the next surface the navigator intends to reach
  template <template <typename> class surface_t = tracking_surface>
  DETRAY_HOST_DEVICE constexpr auto next_surface() const {
    return surface_t<detector_t>{*m_detector, cast_impl().target().surface()};
  }

  /// @param volume_t the volume interface type that is required
  /// @returns the next volume the navigator intends to reach
  template <template <typename> class volume_t = tracking_volume>
  DETRAY_HOST_DEVICE constexpr auto next_volume() const {
    return volume_t<detector_t>{*m_detector,
                                cast_impl().target().surface().volume()};
  }

  /// Navigation reaches final target or leaves detector world. Stop
  /// navigation.
  DETRAY_HOST_DEVICE constexpr void exit() {
    m_status = navigation::status::e_exit;
    DETRAY_VERBOSE_HOST_DEVICE("Exited");
    cast_impl().run_inspector({}, point3_t{0.f, 0.f, 0.f},
                              vector3_t{0.f, 0.f, 0.f}, "Exited: ");
  }

  /// Navigation is being paused by actor: Maintain the navigation state
  /// and resume later
  DETRAY_HOST_DEVICE constexpr void pause() {
    DETRAY_VERBOSE_HOST_DEVICE("Paused by actor");
    cast_impl().run_inspector({}, point3_t{0.f, 0.f, 0.f},
                              vector3_t{0.f, 0.f, 0.f}, "Paused by actor: ");
  }

  /// Navigation state that cannot be recovered from. Leave the other
  /// data for inspection.
  ///
  /// @param custom_msg additional information on the reason for the error
  DETRAY_HOST_DEVICE constexpr void abort(
      const char *custom_msg = "Navigator (unknown reason)") {
    m_status = navigation::status::e_abort;

    /// Wrapper around the custom message that a print inspector can
    /// understand
    struct message_wrapper {
      const char *const m_msg{nullptr};

      DETRAY_HOST_DEVICE
      constexpr const char *operator()() const { return m_msg; }
    };

    assert(custom_msg != nullptr);
    DETRAY_ERROR_HOST("Aborted: " << custom_msg);

    cast_impl().run_inspector({}, point3_t{0.f, 0.f, 0.f},
                              vector3_t{0.f, 0.f, 0.f},
                              "Aborted: ", message_wrapper{custom_msg});
  }

  /// Navigation state that cannot be recovered from. Leave the other
  /// data for inspection.
  ///
  /// @param debug_msg_generator functor that returns additional
  ///                            information on the reason for the error
  template <typename debug_msg_generator_t>
    requires(!std::same_as<char *, debug_msg_generator_t>)
  DETRAY_HOST_DEVICE constexpr void abort(
      const debug_msg_generator_t &debug_msg_generator) {
    m_status = navigation::status::e_abort;

    DETRAY_ERROR_HOST("Aborted: " << debug_msg_generator());

    cast_impl().run_inspector({}, point3_t{0.f, 0.f, 0.f},
                              vector3_t{0.f, 0.f, 0.f},
                              "Aborted: ", debug_msg_generator);
  }

  /// @returns the navigation inspector - const
  DETRAY_HOST_DEVICE
  constexpr const auto &inspector() const { return m_inspector; }

  /// @returns the navigation inspector
  DETRAY_HOST_DEVICE
  constexpr auto &inspector() { return m_inspector; }

  template <typename callable_t>
  DETRAY_HOST_DEVICE constexpr void for_each_valid(callable_t &&f) const {
    for (std::size_t i = this->valid_begin(); i < this->valid_end(); ++i) {
      f(this->candidate_at(i));
    }
  }

 protected:
  /// @returns the first valid cache slot: the current candidate is only
  /// valid while the track is on its surface
  DETRAY_HOST_DEVICE
  constexpr std::size_t valid_begin() const {
    return (is_on_surface() && (m_first == 1)) ? 0u : m_first;
  }

  /// @returns the end of the valid cache range
  DETRAY_HOST_DEVICE
  constexpr std::size_t valid_end() const { return m_n_valid; }

  /// Set the number of valid cache slots, counted from slot 0
  DETRAY_HOST_DEVICE
  constexpr void n_valid(dindex n) {
    m_n_valid = n;
    assert(m_n_valid <= static_cast<dindex>(usable_capacity()));
  }

  /// @returns the number of cache slots that can currently be used. Every
  /// visited candidate beyond the current one is dropped from the cache but
  /// keeps its slot reserved, so the usable depth shrinks by one per surface
  /// until the next (re-)initialization
  DETRAY_HOST_DEVICE
  constexpr std::size_t usable_capacity() const {
    return k_cache_capacity - static_cast<std::size_t>(m_n_dead);
  }

  /// Set the slot of the target: 0 if there is no current candidate, 1 if
  /// the current candidate occupies slot 0
  DETRAY_HOST_DEVICE
  constexpr void first_index(dist_t pos) {
    assert(pos == 0 || pos == 1);
    m_first = pos;
  }

  /// Make the first valid candidate the target (may be the current one)
  DETRAY_HOST_DEVICE
  constexpr void set_next_to_begin() {
    first_index(static_cast<dist_t>(valid_begin()));
    assert(m_first <= static_cast<dist_t>(m_n_valid));
  }

  /// Helper to evict all unreachable/invalid candidates from the cache:
  DETRAY_HOST_DEVICE
  constexpr void evict_invalid() {
    // Depends on previous invalidation of unreachable candidates!
    auto not_reachable = [](const intersection_t &candidate) {
      return candidate.path() == std::numeric_limits<scalar_t>::max();
    };

    dindex i = 0;

    for (; i < k_cache_capacity; ++i) {
      if (not_reachable(this->candidate_at(i))) {
        break;
      }
    }

    n_valid(i);
  }

  template <typename callable_t>
  DETRAY_HOST_DEVICE constexpr void for_each_valid(callable_t &&f) {
    for (std::size_t i = this->valid_begin(); i < this->valid_end(); ++i) {
      auto cand = this->candidate_at(i);
      f(cand);
      this->set_candidate_at(i, cand);
    }
  }

  /// The target has been reached: it becomes the current candidate (slot 0)
  /// and the candidate behind it the new target (slot 1)
  DETRAY_HOST_DEVICE
  constexpr auto advance() -> void {
    assert(n_candidates() >= 1u);
    if (m_first == 0) {
      // The target is already in slot 0
      m_first = 1;
    } else {
      // Drop the previous current candidate and move everything down
      for (std::size_t i = 0u; i + 1u < k_cache_capacity; ++i) {
        set_candidate_at(i, candidate_at(i + 1u));
      }
      auto last = candidate_at(k_cache_capacity - 1u);
      last.set_path(std::numeric_limits<scalar_t>::max());
      set_candidate_at(k_cache_capacity - 1u, last);
      --m_n_valid;
      // The dropped candidate keeps its slot reserved
      ++m_n_dead;
      assert(static_cast<std::size_t>(m_n_dead) < k_cache_capacity);
    }
  }

  DETRAY_HOST_DEVICE
  constexpr void set_current(const candidate_t &c) {
    assert(m_first == 1);
    m_candidates[0] = c;
  };

  DETRAY_HOST_DEVICE
  constexpr void set_target(const candidate_t &c) {
    if (m_first == 0) {
      m_candidates[0] = c;
    } else {
      m_candidates[1] = c;
    }
  };

  DETRAY_HOST_DEVICE
  constexpr void set_candidate_at(std::size_t i, const candidate_t &c) {
    m_candidates[i] = c;
  };

  DETRAY_HOST_DEVICE
  constexpr auto candidate_at(std::size_t i) const -> const candidate_t {
    return m_candidates[i];
  }

  /// Set the status to @param s
  DETRAY_HOST_DEVICE
  constexpr void status(navigation::status s) {
    DETRAY_DEBUG_HOST("Setting nav. status: " << s);
    m_status = s;
  }

  /// Reset the trustlevel to @param t
  DETRAY_HOST_DEVICE
  constexpr void trust_level(navigation::trust_level t) {
    DETRAY_DEBUG_HOST("Setting trust level: " << t);
    m_trust_level = t;
  }

  /// Clear the state
  DETRAY_HOST_DEVICE
  constexpr void clear_cache() {
    // Mark all data in the cache as unreachable
    for (std::size_t i = 0u; i < k_cache_capacity; ++i) {
      auto cand = this->candidate_at(i);
      cand.set_path(std::numeric_limits<scalar_t>::max());
      this->set_candidate_at(i, cand);
    }
    m_first = 0;
    m_n_valid = 0;
    m_n_dead = 0;
  }

  /// Call the navigation inspector
  DETRAY_HOST_DEVICE constexpr void run_inspector(
      [[maybe_unused]] const navigation::config &cfg,
      [[maybe_unused]] const point3_t &track_pos,
      [[maybe_unused]] const vector3_t &track_dir,
      [[maybe_unused]] const char *message) {
    [[maybe_unused]] auto derived = static_cast<const derived_t &>(*this);

    if constexpr (!std::is_same_v<inspector_t, navigation::void_inspector>) {
      cast_impl().inspector()(derived, cfg, track_pos, track_dir, message);
    }

    DETRAY_DEBUG_HOST("" << message << "\n"
                         << detray::navigation::print_state(derived)
                         << detray::navigation::print_candidates(
                                derived, cfg, track_pos, track_dir));
  }

  /// Call the navigation inspector
  template <typename debug_msg_generator_t>
  DETRAY_HOST_DEVICE constexpr void run_inspector(
      [[maybe_unused]] const navigation::config &cfg,
      [[maybe_unused]] const point3_t &track_pos,
      [[maybe_unused]] const vector3_t &track_dir,
      [[maybe_unused]] const char *message,
      [[maybe_unused]] const debug_msg_generator_t &msg_gen) {
    [[maybe_unused]] auto derived = static_cast<const derived_t &>(*this);

    if constexpr (!std::is_same_v<inspector_t, navigation::void_inspector>) {
      cast_impl().inspector()(derived, cfg, track_pos, track_dir, message,
                              msg_gen);
    }

    DETRAY_DEBUG_HOST("" << message << msg_gen() << "\n"
                         << detray::navigation::print_state(derived)
                         << detray::navigation::print_candidates(
                                derived, cfg, track_pos, track_dir));
  }

 private:
  /// @returns a string stream that prints the navigation state details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const derived_t &s) {
    os << detray::navigation::print_state(s)
       << detray::navigation::print_candidates(s, {}, point3_t{0.f, 0.f, 0.f},
                                               vector3_t{0.f, 0.f, 0.f});
    return os;
  }

  /// Our cache of candidates (intersections with any kind of surface).
  /// Slot 0 holds the current candidate once one exists (m_first == 1),
  /// the target sits in slot m_first, further candidates follow in order
  candidate_cache_t m_candidates{};

  /// Detector pointer
  const detector_t *m_detector{nullptr};

  /// The navigation status
  navigation::status m_status{navigation::status::e_unknown};

  /// The navigation direction
  navigation::direction m_direction{navigation::direction::e_forward};

  /// The navigation trust level determines how this states cache is to
  /// be updated in the current navigation call
  navigation::trust_level m_trust_level{navigation::trust_level::e_no_trust};

  /// External mask tolerance, that models noise during track transport
  scalar_t m_external_mask_tol{0.f * unit<scalar_t>::mm};

  /// The number of valid cache slots, counted from slot 0 (includes the
  /// current candidate, if there is one): m_first <= m_n_valid
  dindex m_n_valid{0};

  /// Slot of the target: 0 while there is no current candidate, 1 afterwards
  dist_t m_first{0};

  /// Number of visited candidates that were dropped from the cache since the
  /// last (re-)initialization; their slots stay reserved (usable_capacity())
  dist_t m_n_dead{0};

  /// Index in the detector volume container of current navigation volume
  nav_link_t m_volume_index{0u};

  /// The inspector type of this navigation engine
  DETRAY_NO_UNIQUE_ADDRESS inspector_t m_inspector;
};

}  // namespace detray::navigation
