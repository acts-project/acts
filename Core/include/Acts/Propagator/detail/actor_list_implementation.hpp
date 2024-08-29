// This file is part of the Acts project.
//
// Copyright (C) 2018-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/ActorConcepts.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

#include <utility>

namespace Acts::detail {

namespace {

struct actor_caller {
  template <typename actor_t, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void act(const actor_t& actor, propagator_state_t& state,
                  const stepper_t& stepper, const navigator_t& navigator,
                  Args&&... args) {
    static_assert(
        Actor<actor_t, propagator_state_t, stepper_t, navigator_t, Args...>,
        "The actor does not fulfill the Actor concept");

    if constexpr (ActorHasAct<actor_t, propagator_state_t, stepper_t,
                              navigator_t, Args...>) {
      if constexpr (ActorHasResult<actor_t>) {
        actor.act(state, stepper, navigator,
                  state.template get<detail::result_type_t<actor_t>>(),
                  std::forward<Args>(args)...);
      } else {
        actor.act(state, stepper, navigator, std::forward<Args>(args)...);
      }
    }
  }

  template <typename actor_t, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const actor_t& actor, propagator_state_t& state,
                    const stepper_t& stepper, const navigator_t& navigator,
                    Args&&... args) {
    static_assert(
        Actor<actor_t, propagator_state_t, stepper_t, navigator_t, Args...>,
        "The actor does not fulfill the Actor concept");

    if constexpr (ActorHasAbort<actor_t, propagator_state_t, stepper_t,
                                navigator_t, Args...>) {
      if constexpr (ActorHasResult<actor_t>) {
        return actor.check(state, stepper, navigator,
                           state.template get<detail::result_type_t<actor_t>>(),
                           std::forward<Args>(args)...);
      } else {
        return actor.check(state, stepper, navigator,
                           std::forward<Args>(args)...);
      }
    }

    return false;
  }
};

}  // end of anonymous namespace

/// The dummy list call implementation
template <typename... actors>
struct actor_list_impl;

/// The action list call implementation
/// - it calls 'action' on the current entry of the tuple
/// - then broadcasts the action call to the remaining tuple
template <typename first, typename... others>
struct actor_list_impl<first, others...> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void act(const T& actor_tuple, propagator_state_t& state,
                  const stepper_t& stepper, const navigator_t& navigator,
                  Args&&... args) {
    const auto& this_actor = std::get<first>(actor_tuple);
    actor_caller::act(this_actor, state, stepper, navigator, args...);
    actor_list_impl<others...>::act(actor_tuple, state, stepper, navigator,
                                    args...);
  }

  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const T& actor_tuple, propagator_state_t& state,
                    const stepper_t& stepper, const navigator_t& navigator,
                    Args&&... args) {
    const auto& this_actor = std::get<first>(actor_tuple);
    return actor_caller::check(this_actor, state, stepper, navigator,
                               args...) ||
           actor_list_impl<others...>::check(actor_tuple, state, stepper,
                                             navigator, args...);
  }
};

/// The action list call implementation
/// - it calls 'action' on the last entry of the tuple
template <typename last>
struct actor_list_impl<last> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void act(const T& actors_tuple, propagator_state_t& state,
                  const stepper_t& stepper, const navigator_t& navigator,
                  Args&&... args) {
    const auto& this_actor = std::get<last>(actors_tuple);
    actor_caller::act(this_actor, state, stepper, navigator,
                      std::forward<Args>(args)...);
  }

  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const T& actors_tuple, propagator_state_t& state,
                    const stepper_t& stepper, const navigator_t& navigator,
                    Args&&... args) {
    const auto& this_actor = std::get<last>(actors_tuple);
    return actor_caller::check(this_actor, state, stepper, navigator,
                               std::forward<Args>(args)...);
  }
};

/// The empty action list call implementation
template <>
struct actor_list_impl<> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void act(const T& /*actors_tuple*/, propagator_state_t& /*state*/,
                  const stepper_t& /*stepper*/,
                  const navigator_t& /*navigator*/, Args&&... /*args*/) {}

  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const T& /*actors_tuple*/, propagator_state_t& /*state*/,
                    const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, Args&&... /*args*/) {
    return false;
  }
};

}  // namespace Acts::detail
