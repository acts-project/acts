// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/ActorConcepts.hpp"
#include "Acts/Utilities/Result.hpp"

#include <utility>

namespace Acts::detail {

/// Call the act method of an actor, with its state if it has one
///
/// An actor without an act method does nothing.
template <typename actor_t, typename actor_state_t, typename propagator_state_t,
          typename stepper_t, typename navigator_t, typename... Args>
Result<void> callActorAct(const actor_t& actor, actor_state_t&& actorState,
                          propagator_state_t& state, const stepper_t& stepper,
                          const navigator_t& navigator, Args&&... args)
  requires(Actor<actor_t, propagator_state_t, stepper_t, navigator_t, Args...>)
{
  if constexpr (ActorHasActWithoutState<actor_t, propagator_state_t, stepper_t,
                                        navigator_t, Args...>) {
    static_cast<void>(actorState);
    return actor.act(state, stepper, navigator, std::forward<Args>(args)...);
  } else if constexpr (ActorHasActWithState<actor_t, propagator_state_t,
                                            stepper_t, navigator_t, Args...>) {
    return actor.act(state, stepper, navigator, actorState,
                     std::forward<Args>(args)...);
  } else {
    static_cast<void>(actorState);
    return Result<void>::success();
  }
}

/// Call the checkAbort method of an actor, with its state if it has one
///
/// An actor without a checkAbort method never aborts.
template <typename actor_t, typename actor_state_t, typename propagator_state_t,
          typename stepper_t, typename navigator_t, typename... Args>
bool callActorCheckAbort(const actor_t& actor, actor_state_t&& actorState,
                         propagator_state_t& state, const stepper_t& stepper,
                         const navigator_t& navigator, Args&&... args)
  requires(Actor<actor_t, propagator_state_t, stepper_t, navigator_t, Args...>)
{
  if constexpr (ActorHasAbortWithoutState<actor_t, propagator_state_t,
                                          stepper_t, navigator_t, Args...>) {
    static_cast<void>(actorState);
    return actor.checkAbort(state, stepper, navigator,
                            std::forward<Args>(args)...);
  } else if constexpr (ActorHasAbortWithState<actor_t, propagator_state_t,
                                              stepper_t, navigator_t,
                                              Args...>) {
    return actor.checkAbort(state, stepper, navigator, actorState,
                            std::forward<Args>(args)...);
  } else {
    static_cast<void>(actorState);
    return false;
  }
}

}  // namespace Acts::detail
