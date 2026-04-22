// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/composite_actor.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <concepts>
#include <optional>
#include <type_traits>
#include <utility>

namespace detray {

/// The interface to the actors and aborters in the propagation.
///
/// It can hold both simple actors, as well as an actor with its observers.
/// The states of the actors need to be passed to the chain in an external tuple
///
/// @tparam actors_t the types of the actors in the chain.
template <concepts::actor... actors_t>
class actor_chain {
 public:
  /// Types of the actors that are registered in the chain
  using actor_tuple = dtuple<actors_t...>;

  // Tuple of actor states (including states of observing actors, if present)
  using state_tuple =
      detail::unique_t<detail::tuple_cat_t<detail::state_tuple_t<actors_t>...>>;

  // Tuple of state references that is used in the propagator
  using state_ref_tuple = detail::unique_t<
      detail::tuple_cat_t<detail::state_ref_tuple_t<actors_t>...>>;

  /// Call all actors in the chain.
  ///
  /// @param states the states of the actors.
  /// @param p_state the propagation state.
  template <typename actor_states_t, typename propagator_state_t>
  DETRAY_HOST_DEVICE void operator()(actor_states_t &states,
                                     propagator_state_t &p_state) const {
    DETRAY_VERBOSE_HOST_DEVICE("Run actors...");
    run(states, p_state, std::make_index_sequence<sizeof...(actors_t)>{});
  }

  /// @returns the actor list
  DETRAY_HOST_DEVICE constexpr const actor_tuple &actors() const {
    return m_actors;
  }

  /// @returns a tuple of default constructible actor states
  DETRAY_HOST_DEVICE
  static constexpr auto make_default_actor_states() {
    // Only possible if each state is default initializable (including
    // obsevers to the actors in actors_t)
    if constexpr (detail::tuple_all_v<std::is_default_constructible,
                                      state_tuple>) {
      return state_tuple();
    } else {
      return std::nullopt;
    }
  }

  /// @returns a tuple of reference for every state in the tuple @param t
  DETRAY_HOST_DEVICE static constexpr state_ref_tuple setup_actor_states(
      state_tuple &t) {
    return setup_actor_states(
        t, std::make_index_sequence<detail::tuple_size_v<state_tuple>>{});
  }

 private:
  /// Call the actors. Either single actor or composition.
  ///
  /// @param actr the actor (might be a composite actor)
  /// @param states states of all actors (only bare actors)
  /// @param p_state the state of the propagator (stepper and navigator)
  template <concepts::actor actor_t, typename actor_states_t,
            typename propagator_state_t>
  DETRAY_HOST_DEVICE inline void run(const actor_t &actr,
                                     actor_states_t &states,
                                     propagator_state_t &p_state) const {
    if constexpr (!typename actor_t::is_comp_actor()) {
      // No actor state defined (empty)
      if constexpr (std::same_as<typename actor_t::state,
                                 detray::base_actor::state>) {
        actr(p_state);
      } else {
        actr(detail::get<typename actor_t::state &>(states), p_state);
      }
    } else {
      actr(states, p_state);
    }
  }

  /// Resolve the actor calls.
  ///
  /// @param states states of all actors
  /// @param p_state the state of the propagator (stepper and navigator)
  template <typename actor_states_t, typename propagator_state_t,
            std::size_t... indices>
  DETRAY_HOST_DEVICE inline void run(
      actor_states_t &states, propagator_state_t &p_state,
      std::index_sequence<indices...> /*ids*/) const {
    (run(detail::get<indices>(m_actors), states, p_state), ...);
  }

  /// @returns a tuple of reference for every state in the tuple @param t
  template <std::size_t... indices>
  DETRAY_HOST_DEVICE static constexpr state_ref_tuple setup_actor_states(
      state_tuple &t, std::index_sequence<indices...> /*ids*/) {
    return detray::tie(detail::get<indices>(t)...);
  }

  /// Tuple of actors
  DETRAY_NO_UNIQUE_ADDRESS actor_tuple m_actors = {};
};

/// Empty actor chain (placeholder)
template <>
class actor_chain<> {
 public:
  using actor_tuple = dtuple<>;
  using state_tuple = dtuple<>;
  using state_ref_tuple = dtuple<>;

  /// Empty states replaces a real actor states container
  using state = state_tuple;

  /// Call to actors does nothing.
  ///
  /// @param states the states of the actors.
  /// @param p_state the propagation state.
  template <typename actor_states_t, typename propagator_state_t>
  DETRAY_HOST_DEVICE constexpr void operator()(
      actor_states_t & /*states*/, propagator_state_t & /*p_state*/) const {
    /*Do nothing*/
  }

  /// @returns the actor list
  DETRAY_HOST_DEVICE constexpr const actor_tuple &actors() const {
    return m_actors;
  }

  /// @returns an empty state
  DETRAY_HOST_DEVICE
  static consteval state_tuple make_default_actor_states() {
    return dtuple<>{};
  }

  /// @returns an empty state
  DETRAY_HOST_DEVICE static constexpr state_ref_tuple setup_actor_states(
      const state_tuple &) {
    return {};
  }

 private:
  DETRAY_NO_UNIQUE_ADDRESS actor_tuple m_actors = {};
};

}  // namespace detray
