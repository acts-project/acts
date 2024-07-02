// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/detail/MPL/type_collector.hpp"

#include <utility>

namespace Acts::detail {

namespace {

/// The action caller struct, it's called with the the right result object,
/// which is taken out from the result tuple
template <bool has_result = true>
struct action_caller {
  template <typename actor_t, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void action(const actor_t& act, propagator_state_t& state,
                     const stepper_t& stepper, const navigator_t& navigator,
                     Args&&... args) {
    act(state, stepper, navigator,
        state.template get<detail::result_type_t<actor_t>>(),
        std::forward<Args>(args)...);
  }
};

/// The action caller struct, without result object
template <>
struct action_caller<false> {
  template <typename actor_t, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void action(const actor_t& act, propagator_state_t& state,
                     const stepper_t& stepper, const navigator_t& navigator,
                     Args&&... args) {
    act(state, stepper, navigator, std::forward<Args>(args)...);
  }
};
}  // end of anonymous namespace

/// The dummy list call implementation
template <typename... actors>
struct action_list_impl;

/// The action list call implementation
/// - it calls 'action' on the current entry of the tuple
/// - then broadcasts the action call to the remaining tuple
template <typename first, typename... others>
struct action_list_impl<first, others...> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void action(const T& actors_tuple, propagator_state_t& state,
                     const stepper_t& stepper, const navigator_t& navigator,
                     Args&&... args) {
    constexpr bool has_result = has_result_type_v<first>;
    const auto& this_actor = std::get<first>(actors_tuple);
    action_caller<has_result>::action(this_actor, state, stepper, navigator,
                                      args...);
    action_list_impl<others...>::action(actors_tuple, state, stepper, navigator,
                                        args...);
  }
};

/// The action list call implementation
/// - it calls 'action' on the last entry of the tuple
template <typename last>
struct action_list_impl<last> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void action(const T& actors_tuple, propagator_state_t& state,
                     const stepper_t& stepper, const navigator_t& navigator,
                     Args&&... args) {
    constexpr bool has_result = has_result_type_v<last>;
    const auto& this_actor = std::get<last>(actors_tuple);
    action_caller<has_result>::action(this_actor, state, stepper, navigator,
                                      std::forward<Args>(args)...);
  }
};

/// The empty action list call implementation
template <>
struct action_list_impl<> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static void action(const T& /*actors_tuple*/, propagator_state_t& /*state*/,
                     const stepper_t& /*stepper*/,
                     const navigator_t& /*navigator*/, Args&&... /*args*/) {}
};

}  // namespace Acts::detail
