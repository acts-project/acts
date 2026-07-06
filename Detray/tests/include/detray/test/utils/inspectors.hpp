// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: Remove this when gcc fixes their false positives.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#endif

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/math.hpp"
#include "detray/navigation/detail/print_state.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/propagator/detail/print_stepper_state.hpp"
#include "detray/propagator/stepping_config.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace detray {

/// An inspector that aggregates a number of different inspectors.
template <typename... Inspectors>
struct aggregate_inspector {
  using view_type = dmulti_view<detail::get_view_t<Inspectors>...>;
  using const_view_type = dmulti_view<detail::get_view_t<const Inspectors>...>;

  using inspector_tuple_t = std::tuple<Inspectors...>;
  inspector_tuple_t _inspectors{};

  /// Default constructor
  aggregate_inspector() = default;

  /// Construct from the inspector @param view type. Mainly used device-side.
  template <concepts::device_view view_t>
  DETRAY_HOST_DEVICE explicit aggregate_inspector(view_t &view)
      : _inspectors(unroll_views(
            view, std::make_index_sequence<sizeof...(Inspectors)>{})) {}

  /// Inspector interface
  template <unsigned int current_id = 0, typename state_type,
            concepts::point3D point3_t, concepts::vector3D vector3_t,
            typename... Args>
  DETRAY_HOST_DEVICE auto operator()(state_type &state,
                                     const navigation::config &cfg,
                                     const point3_t &pos, const vector3_t &dir,
                                     const char *message, Args &&...args) {
    // Call inspector
    std::get<current_id>(_inspectors)(state, cfg, pos, dir, message,
                                      std::forward<Args>(args)...);

    // Next inspector
    if constexpr (current_id < std::tuple_size<inspector_tuple_t>::value - 1) {
      return operator()<current_id + 1>(state, cfg, pos, dir, message,
                                        std::forward<Args>(args)...);
    }
  }

  /// Inspector interface
  template <unsigned int current_id = 0, typename state_type, typename... Args>
  DETRAY_HOST_DEVICE auto operator()(state_type &state, const char *message,
                                     Args &&...args) {
    // Call inspector
    std::get<current_id>(_inspectors)(state, message);

    // Next inspector
    if constexpr (current_id < std::tuple_size<inspector_tuple_t>::value - 1) {
      return operator()<current_id + 1>(state, message,
                                        std::forward<Args>(args)...);
    }
  }

  /// @returns a specific inspector by type
  template <typename inspector_t>
  DETRAY_HOST_DEVICE constexpr decltype(auto) get() {
    return std::get<inspector_t>(_inspectors);
  }

  /// @returns a specific inspector by type - const
  template <typename inspector_t>
  DETRAY_HOST_DEVICE constexpr decltype(auto) get() const {
    return std::get<inspector_t>(_inspectors);
  }

  /// @returns a specific inspector by index
  template <std::size_t I>
  DETRAY_HOST_DEVICE constexpr decltype(auto) get() {
    return std::get<I>(_inspectors);
  }

  /// @returns a specific inspector by index - const
  template <std::size_t I>
  DETRAY_HOST_DEVICE constexpr decltype(auto) get() const {
    return std::get<I>(_inspectors);
  }

  /// @returns a tuple constructed from the inspector @param view s.
  template <concepts::device_view view_t, std::size_t... I>
  DETRAY_HOST_DEVICE constexpr auto unroll_views(
      view_t &view, std::index_sequence<I...> /*seq*/) {
    return detail::make_tuple<std::tuple>(
        Inspectors(detail::get<I>(view.m_view))...);
  }
};

namespace navigation {

namespace detail {

/// Record of a surface intersection along a track
template <typename intersetion_t>
struct candidate_record {
  using algebra_type = typename intersetion_t::algebra_type;
  using scalar_type = dscalar<algebra_type>;
  using point3_type = dpoint3D<algebra_type>;
  using vector3_type = dvector3D<algebra_type>;
  using intersection_type = intersetion_t;

  constexpr candidate_record() = default;

  /// The particle charge is not known in the navigation, but might be
  /// provided in a different context
  DETRAY_HOST_DEVICE
  constexpr candidate_record(
      const point3_type &position, const vector3_type &direction,
      const intersection_type &intr,
      const scalar_type q = detray::detail::invalid_value<scalar_type>(),
      const scalar_type p = detray::detail::invalid_value<scalar_type>())
      : pos{position},
        dir{direction},
        intersection{intr},
        charge{q},
        p_mag{p} {}

  /// Current global track position
  point3_type pos{0.f, 0.f, 0.f};
  /// Current global track direction
  vector3_type dir{0.f, 0.f, 1.f};
  /// The intersection result
  intersetion_t intersection{};
  /// Charge hypothesis of the particle (invalid value if not known)
  scalar_type charge{detray::detail::invalid_value<scalar_type>()};
  /// Current momentum magnitude of the particle
  scalar_type p_mag{1.f};
};

}  // namespace detail

/// A navigation inspector that relays information about the encountered
/// objects whenever the navigator reaches one or more status flags
template <typename candidate_t, template <typename...> class vector_t = dvector,
          status... navigation_status>
struct object_tracer {
  using candidate_record_t = detail::candidate_record<candidate_t>;
  using scalar_t = typename candidate_record_t::scalar_type;

  using view_type = dvector_view<candidate_record_t>;
  using const_view_type = dvector_view<const candidate_record_t>;

  /// Default constructor
  object_tracer() = default;

  /// Device-side construction from a vecmem based view type
  DETRAY_HOST_DEVICE explicit object_tracer(
      dvector_view<candidate_record_t> &view)
      : object_trace(view) {}

  // Record all objects the navigator encounters
  vector_t<candidate_record_t> object_trace;
  dindex current_vol{dindex_invalid};
  const scalar_t inv_pos{detray::detail::invalid_value<scalar_t>()};
  typename candidate_record_t::point3_type last_pos = {inv_pos, inv_pos,
                                                       inv_pos};
  typename candidate_record_t::vector3_type last_dir = {0.f, 0.f, 0.f};

  /// Inspector interface
  template <typename state_type, concepts::point3D point3_t,
            concepts::vector3D vector3_t, typename... Args>
  DETRAY_HOST_DEVICE auto operator()(const state_type &state,
                                     const navigation::config & /*unused*/,
                                     const point3_t &pos, const vector3_t &dir,
                                     const char * /*message*/,
                                     Args &&.../*unused*/) {
    DETRAY_VERBOSE_HOST_DEVICE("Actor: Check navigation status...");

    // Record the candidate of an encountered object
    if ((is_status(state.status(), navigation_status) || ...)) {
      DETRAY_VERBOSE_HOST_DEVICE("Actor: ...Status matched:");
      // Reached a new position: log it
      // Also log volume switches that happen without position update
      const bool first_obj{(last_pos == point3_t{inv_pos, inv_pos, inv_pos}) &&
                           (last_dir == point3_t{0.f, 0.f, 0.f})};

      if ((state.is_on_portal() && current_vol != state.volume()) ||
          object_trace.empty() || first_obj ||
          state.current().surface().identifier() !=
              object_trace.back().intersection.surface().identifier()) {
        DETRAY_VERBOSE_HOST_DEVICE("Actor: -> Record surface %d",
                                   state.current_surface().index());

        object_trace.push_back({pos, dir, state.current()});
        last_pos = pos;
        last_dir = dir;
        current_vol = state.volume();
      }
    }
  }

  /// Inspector interface
  template <typename state_type>
  DETRAY_HOST_DEVICE auto operator()(const state_type & /*state*/,
                                     const char * /*message*/) { /* Do nothing*/
  }

  /// @returns a specific candidate from the trace
  DETRAY_HOST_DEVICE
  constexpr const candidate_record_t &operator[](std::size_t i) const {
    return object_trace[i];
  }

  /// @returns the full object trace
  DETRAY_HOST_DEVICE
  constexpr const auto &trace() const { return object_trace; }

  /// Compares a navigation status with the tracers references
  DETRAY_HOST_DEVICE
  constexpr bool is_status(const status &nav_stat, const status &ref_stat) {
    return (nav_stat == ref_stat);
  }
};

/// A navigation inspector that prints information about the current navigation
/// state. Meant for debugging.
struct print_inspector {
  using view_type = dvector_view<char>;
  using const_view_type = dvector_view<const char>;

  struct void_generator {};

  /// Default constructor
  print_inspector() = default;

  /// Copy constructor ensures that the string stream is set up identically
  print_inspector(const print_inspector &other)
      : debug_stream(other.debug_stream.str()) {}

  /// Move constructor
  print_inspector(print_inspector &&other) = default;

  /// Default destructor
  ~print_inspector() = default;

  /// Copy assignemten operator ensures that the string stream is set up
  /// identically
  print_inspector &operator=(const print_inspector &other) {
    // Reset
    debug_stream.str(std::string());
    debug_stream.clear();
    // Move new debug string in
    debug_stream << other.debug_stream.str();
    return *this;
  }

  /// Move assignemten operator
  print_inspector &operator=(print_inspector &&other) = default;

  /// Inspector interface. Gathers detailed information during navigation
  template <typename state_type, concepts::point3D point3_t,
            concepts::vector3D vector3_t,
            typename message_generator_t = void_generator>
  auto operator()(const state_type &state, const navigation::config &cfg,
                  const point3_t &track_pos, const vector3_t &track_dir,
                  const char *message,
                  const message_generator_t &msg_gen = {}) {
    std::string msg(message);
    debug_stream << msg;
    if constexpr (!std::same_as<message_generator_t, void_generator>) {
      debug_stream << msg_gen();

      if (state.status() == navigation::status::e_abort) {
        fata_error_msg = msg_gen();
      }
    }
    debug_stream << std::endl;
    debug_stream << "----------------------------------------" << std::endl;

    debug_stream << navigation::print_state(state);
    debug_stream << navigation::print_candidates(state, cfg, track_pos,
                                                 track_dir);

    debug_stream << std::endl << std::endl;
  }

  /// Inspector interface. Print basic state information
  template <typename state_type, typename message_generator_t = void_generator>
  auto operator()(const state_type &state, const char *message,
                  const message_generator_t &msg_gen = {}) {
    std::string msg(message);
    debug_stream << msg;
    if constexpr (!std::same_as<message_generator_t, void_generator>) {
      debug_stream << msg_gen();

      if (state.status() == navigation::status::e_abort) {
        fata_error_msg = msg_gen();
      }
    }
    debug_stream << std::endl;
    debug_stream << "----------------------------------------" << std::endl;

    debug_stream << navigation::print_state(state);

    debug_stream << std::endl << std::endl;
  }

  /// @returns a string representation of the gathered information
  std::string to_string() const { return debug_stream.str(); }

  /// Gathers navigation information across navigator update calls
  std::stringstream debug_stream{};
  /// Special message that is collected if the navigator hits a fatal error
  std::string fata_error_msg{""};
};

}  // namespace navigation

namespace stepping {

/// A stepper inspector that prints information about the current stepper
/// state. Meant for debugging.
struct print_inspector {
  /// Default constructor
  print_inspector() = default;

  /// Copy constructor ensures that the string stream is set up identically
  print_inspector(const print_inspector &other)
      : debug_stream(other.debug_stream.str()) {}

  /// Copy assignemten operator ensures that the string stream is set up
  /// identically
  print_inspector &operator=(const print_inspector &other) {
    // Reset
    debug_stream.str(std::string());
    debug_stream.clear();
    // Move new debug string in
    debug_stream << other.debug_stream.str();
    return *this;
  }

  /// Move assignemten operator
  print_inspector &operator=(print_inspector &&other) = default;

  /// Gathers stepping information from inside the stepper methods
  std::stringstream debug_stream{};

  /// Inspector interface. Gathers detailed information during stepping
  template <typename state_type, concepts::scalar scalar_t>
  void operator()(const state_type &state, const stepping::config & /*unused*/,
                  const char *message, const scalar_t dist) {
    std::string msg(message);

    debug_stream << msg << std::endl;
    debug_stream << detray::stepping::print_state(state, dist);
  }

  /// Inspector interface. Gathers detailed information during stepping
  template <typename state_type, concepts::scalar scalar_t>
  void operator()(const state_type &state, const stepping::config & /*unused*/,
                  const char *message, const scalar_t /*dist*/,
                  const std::size_t n_trials, const scalar_t step_scalor) {
    std::string msg(message);

    debug_stream << msg << std::endl;
    debug_stream << detray::stepping::print_state(state, n_trials, step_scalor);
  }

  /// @returns a string representation of the gathered information
  std::string to_string() const { return debug_stream.str(); }
};

}  // namespace stepping

}  // namespace detray
