// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/ActorConcepts.hpp"
#include "Acts/Propagator/detail/ActorCalls.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeList.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"

#include <tuple>
#include <utility>

namespace Acts {

/// @brief Composite actor that runs a list of actors
///
/// The propagator runs a single actor. This actor combines several actors into
/// one. It calls the actors in the order of the list. The act calls stop at
/// the first failure, and the propagation aborts if any actor aborts.
///
/// The state of the composite holds the states of the actors that have one.
/// Access them with `actorState.get<actor_t::state_type>()`.
///
/// @note The composite is fully responsible for the abort conditions of the
///       propagation, like any other actor. Add @c DefaultActor to stop at the
///       target surface and at the path limit. Put it first, so the loop
///       protection sees the start momentum before other actors change it.
template <typename... actors_t>
  requires((std::is_default_constructible_v<actors_t> && ...) &&
           (std::is_copy_constructible_v<actors_t> && ...) &&
           (std::is_move_constructible_v<actors_t> && ...) &&
           ((Types::count<actors_t, TypeList<actors_t...>>::value == 1) && ...))
struct ActorList {
  /// State of the composite: the states of the actors that have one
  using state_type = typename Types::apply<
      detail::Extendable,
      typename Types::map<detail::ActorStateTypeImpl,
                          typename Types::filter<ActorHasStateStruct,
                                                 actors_t...>::type>::type>::
      type;

  /// Const retrieval of an actor of a specific type
  ///
  /// @tparam actor_t Type of the Actor to be retrieved
  /// @return Const reference to the requested actor
  template <typename actor_t>
  const actor_t& get() const {
    return std::get<actor_t>(m_actors);
  }

  /// Non-const retrieval of an actor of a specific type
  ///
  /// @tparam actor_t Type of the Actor to be retrieved
  /// @return Reference to the requested actor
  template <typename actor_t>
  actor_t& get() {
    return std::get<actor_t>(m_actors);
  }

  /// Act call which broadcasts the call to the actors of the list
  ///
  /// @tparam propagator_state_t is the state type of the propagator
  /// @tparam stepper_t Type of the stepper used for the propagation
  /// @tparam navigator_t Type of the navigator used for the propagation
  ///
  /// @param [in,out] state This is the propagator state object
  /// @param [in] stepper The stepper in use
  /// @param [in] navigator The navigator in use
  /// @param [in,out] actorState The state of the composite
  /// @param [in] args The arguments to be passed to the actions
  /// @return A Result<void> indicating success or failure
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                   const navigator_t& navigator, state_type& actorState,
                   Args&&... args) const {
    Result<void> result = Result<void>::success();
    std::apply(
        [&](const actors_t&... actor) {
          // Stop at the first actor that fails
          static_cast<void>(
              ((result = detail::callActorAct(
                    actor, stateOf<actors_t>(actorState), state, stepper,
                    navigator, std::forward<Args>(args)...))
                   .ok() &&
               ...));
        },
        m_actors);
    return result;
  }

  /// Check call which broadcasts the call to the actors of the list
  ///
  /// @tparam propagator_state_t is the state type of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @tparam navigator_t Type of the navigator
  ///
  /// @param [in,out] state is the state object from the propagator
  /// @param [in] stepper Stepper used for the propagation
  /// @param [in] navigator Navigator used for the propagation
  /// @param [in,out] actorState The state of the composite
  /// @param [in] args are the arguments to be passed to the actors
  /// @return True if propagation should be aborted, false otherwise
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  bool checkAbort(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, state_type& actorState,
                  Args&&... args) const {
    return std::apply(
        [&](const actors_t&... actor) {
          return (detail::callActorCheckAbort(
                      actor, stateOf<actors_t>(actorState), state, stepper,
                      navigator, std::forward<Args>(args)...) ||
                  ...);
        },
        m_actors);
  }

 private:
  std::tuple<actors_t...> m_actors;

  /// The state of one actor of the list, or a placeholder if it has none
  template <typename actor_t>
  static decltype(auto) stateOf(state_type& actorState) {
    if constexpr (ActorHasState<actor_t>) {
      return actorState.template get<typename actor_t::state_type>();
    } else {
      return NoActorState{};
    }
  }
};

}  // namespace Acts
