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

// Include all core actors
#include "detray/definitions/pdg_particle.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <concepts>
#include <optional>

namespace detray::concepts {

/// Concept for a simple actor
template <typename A>
concept actor = std::derived_from<A, detray::base_actor> &&
                requires(const A a) { typename A::state; };

/// Concept for an actor including observing actors
template <typename A>
concept composite_actor =
    actor<A> && A::is_comp_actor::value && requires(const A ca) {
      typename A::observer_states;
      typename A::observer_state_refs;
    };

/// Concept for the actor chain that is being run in the propagator
template <typename A>
concept actor_chain = requires(const A c, typename A::state_tuple s) {
  typename A::actor_tuple;
  typename A::state_tuple;
  typename A::state_ref_tuple;

  { c.actors() } -> std::same_as<const typename A::actor_tuple &>;

  {
    c.make_default_actor_states()
  } -> detray::concepts::any_of<typename A::state_tuple, std::nullopt_t>;

  { c.setup_actor_states(s) } -> std::same_as<typename A::state_ref_tuple>;
};

/// Check if a state type belongs to an actor or actor chain
template <typename S, typename T>
concept is_state_of =
    (actor_chain<T> && (detail::is_permutation_v<std::remove_cvref_t<S>,
                                                 typename T::state_ref_tuple> ||
                        detail::is_permutation_v<std::remove_cvref_t<S>,
                                                 typename T::state_tuple>)) ||

    (actor<T> && std::same_as<std::remove_cvref_t<S>, typename T::state>);

/// Check if a type can function as a propagator state
template <typename S>
concept propagator_state = requires {
  typename S::detector_type;
  typename S::detector_type::scalar_type;
  typename S::stepper_state_type;
  typename S::navigator_state_type;
  typename S::context_type;
} && requires(S &s) {
  requires requires(
      const pdg_particle<typename S::detector_type::scalar_type> &ptc) {
    { s.set_particle(ptc) };
  };

  requires requires(bool b) {
    { s.heartbeat(b) };
    { s.debug(b) };
  };

  { s.stepping() } -> std::same_as<typename S::stepper_state_type &>;
  { s.navigation() } -> std::same_as<typename S::navigator_state_type &>;
  { s.context() } -> std::same_as<typename S::context_type &>;
} && requires(const S &s) {
  { s.debug() } -> std::same_as<bool>;
  { s.is_alive() } -> std::same_as<bool>;
  { s.heartbeat() } -> std::same_as<bool>;
  { s.stepping() } -> std::same_as<const typename S::stepper_state_type &>;
  { s.navigation() } -> std::same_as<const typename S::navigator_state_type &>;
  { s.context() } -> std::same_as<const typename S::context_type &>;
};

/// Check if a type can function as a propagator state for a specific
/// propagator type
template <typename S, typename P>
concept is_propagator_state_of = propagator_state<S> && requires {
  typename P::detector_type;
  typename P::stepper_type;
  typename P::navigator_type;

  requires std::same_as<typename P::detector_type, typename S::detector_type>;
  requires std::same_as<typename P::stepper_type::state,
                        typename S::stepper_state_type>;
  requires std::same_as<typename P::navigator_type::state,
                        typename S::navigator_state_type>;
};

/// Check if a type corresponds to the magnetic field type of a stepper
template <typename S, typename F>
concept is_field_of = requires { typename S::magnetic_field_type; } &&
                      std::same_as<F, typename S::magnetic_field_type>;

}  // namespace detray::concepts
