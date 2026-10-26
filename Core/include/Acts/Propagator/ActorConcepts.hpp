// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Result.hpp"

#include <concepts>
#include <type_traits>
#include <utility>

namespace Acts {

template <typename actor_t>
concept ActorHasState = requires { typename actor_t::state_type; };

/// Helper struct exposing ActorHasState as a boolean.
template <typename actor_t>
struct ActorHasStateStruct {
  /// Whether the actor has a state type
  static constexpr bool value = ActorHasState<actor_t>;
};

/// Placeholder for the state of an actor without a state type
struct NoActorState {};

namespace detail {
template <typename actor_t>
struct ActorStateTypeImpl {
  using type = NoActorState;
};
template <ActorHasState actor_t>
struct ActorStateTypeImpl<actor_t> {
  using type = typename actor_t::state_type;
};
}  // namespace detail

/// The state type of an actor, or @c NoActorState if it has none
template <typename actor_t>
using ActorStateType = typename detail::ActorStateTypeImpl<actor_t>::type;

template <typename actor_t, typename propagator_state_t, typename stepper_t,
          typename navigator_t, typename... Args>
concept ActorHasActWithoutState = requires(
    const actor_t& a, propagator_state_t& state, const stepper_t& stepper,
    const navigator_t& navigator, Args&&... args) {
  {
    a.act(state, stepper, navigator, std::forward<Args>(args)...)
  } -> std::same_as<Result<void>>;
};

template <typename actor_t, typename propagator_state_t, typename stepper_t,
          typename navigator_t, typename... Args>
concept ActorHasActWithState =
    ActorHasState<actor_t> &&
    requires(const actor_t& a, propagator_state_t& state,
             const stepper_t& stepper, const navigator_t& navigator,
             typename actor_t::state_type& actorState, Args&&... args) {
      {
        a.act(state, stepper, navigator, actorState,
              std::forward<Args>(args)...)
      } -> std::same_as<Result<void>>;
    };

template <typename actor_t, typename propagator_state_t, typename stepper_t,
          typename navigator_t, typename... Args>
concept ActorHasOldVoidInterface =
    // Check without state parameter (for actors without state_type)
    requires(const actor_t& a, propagator_state_t& state,
             const stepper_t& stepper, const navigator_t& navigator,
             Args&&... args) {
      {
        a.act(state, stepper, navigator, std::forward<Args>(args)...)
      } -> std::same_as<void>;
    } ||
    // Check with state parameter (for actors with state_type)
    (ActorHasState<actor_t> &&
     requires(const actor_t& a, propagator_state_t& state,
              const stepper_t& stepper, const navigator_t& navigator,
              typename actor_t::state_type& actorState, Args&&... args) {
       {
         a.act(state, stepper, navigator, actorState,
               std::forward<Args>(args)...)
       } -> std::same_as<void>;
     });

template <typename actor_t, typename propagator_state_t, typename stepper_t,
          typename navigator_t, typename... Args>
concept ActorHasAct =
    ActorHasActWithoutState<actor_t, propagator_state_t, stepper_t, navigator_t,
                            Args...> ||
    ActorHasActWithState<actor_t, propagator_state_t, stepper_t, navigator_t,
                         Args...>;

template <typename actor_t, typename propagator_state_t, typename stepper_t,
          typename navigator_t, typename... Args>
concept ActorHasAbortWithoutState = requires(
    const actor_t& a, propagator_state_t& state, const stepper_t& stepper,
    const navigator_t& navigator, Args&&... args) {
  {
    a.checkAbort(state, stepper, navigator, std::forward<Args>(args)...)
  } -> std::same_as<bool>;
};

template <typename actor_t, typename propagator_state_t, typename stepper_t,
          typename navigator_t, typename... Args>
concept ActorHasAbortWithState =
    ActorHasState<actor_t> &&
    requires(const actor_t& a, propagator_state_t& state,
             const stepper_t& stepper, const navigator_t& navigator,
             typename actor_t::state_type& actorState, Args&&... args) {
      {
        a.checkAbort(state, stepper, navigator, actorState,
                     std::forward<Args>(args)...)
      } -> std::same_as<bool>;
    };

template <typename actor_t, typename propagator_state_t, typename stepper_t,
          typename navigator_t, typename... Args>
concept ActorHasAbort =
    ActorHasAbortWithoutState<actor_t, propagator_state_t, stepper_t,
                              navigator_t, Args...> ||
    ActorHasAbortWithState<actor_t, propagator_state_t, stepper_t, navigator_t,
                           Args...>;

template <typename actor_t, typename propagator_state_t, typename stepper_t,
          typename navigator_t, typename... Args>
concept Actor = (ActorHasAct<actor_t, propagator_state_t, stepper_t,
                             navigator_t, Args...> ||
                 ActorHasAbort<actor_t, propagator_state_t, stepper_t,
                               navigator_t, Args...>) &&
                !ActorHasOldVoidInterface<actor_t, propagator_state_t,
                                          stepper_t, navigator_t, Args...>;

}  // namespace Acts
