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

/// This is the caller used if the aborter uses the result and that result type
/// exists
template <bool has_result = true>
struct abort_checker {
  template <typename aborter_t, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const aborter_t& aborter, propagator_state_t& state,
                    const stepper_t& stepper, const navigator_t& navigator,
                    Args&&... args) {
    using action_type = action_type_t<aborter_t>;
    using result_type = result_type_t<action_type>;

    return aborter(state, stepper, navigator, state.template get<result_type>(),
                   std::forward<Args>(args)...);
  }
};

/// This is the caller used if the aborter only uses the cache it has no access
/// to the result
template <>
struct abort_checker<false> {
  template <typename aborter_t, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const aborter_t& aborter, propagator_state_t& state,
                    const stepper_t& stepper, const navigator_t& navigator,
                    Args&&... args) {
    return aborter(state, stepper, navigator, std::forward<Args>(args)...);
  }
};
}  // end of anonymous namespace

template <typename... conditions>
struct abort_list_impl;

/// This is the check call on the a list of conditions it calls the apparent
/// condition and broadcasts the call to the remaining ones
template <typename first, typename... others>
struct abort_list_impl<first, others...> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const T& aborters_tuple, propagator_state_t& state,
                    const stepper_t& stepper, const navigator_t& navigator,
                    Args&&... args) {
    // get the right helper for calling the abort condition
    constexpr bool has_result = has_action_type_v<first>;

    // get the cache abort condition
    const auto& this_aborter = std::get<first>(aborters_tuple);
    // - check abort conditions recursively
    // - make use of short-circuit evaluation
    // -> skip remaining conditions if this abort condition evaluates to true

    bool abort = abort_checker<has_result>::check(this_aborter, state, stepper,
                                                  navigator, args...) ||
                 abort_list_impl<others...>::check(aborters_tuple, state,
                                                   stepper, navigator, args...);

    return abort;
  }
};

/// This is the check call on the a last of all conditions
template <typename last>
struct abort_list_impl<last> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const T& aborters_tuple, propagator_state_t& state,
                    const stepper_t& stepper, const navigator_t& navigator,
                    Args&&... args) {
    // get the right helper for calling the abort condition
    constexpr bool has_result = has_action_type_v<last>;
    const auto& this_aborter = std::get<last>(aborters_tuple);
    return abort_checker<has_result>::check(
        this_aborter, state, stepper, navigator, std::forward<Args>(args)...);
  }
};

/// This is the empty call list - never abort
template <>
struct abort_list_impl<> {
  template <typename T, typename propagator_state_t, typename stepper_t,
            typename navigator_t, typename... Args>
  static bool check(const T& /*aborter_tuple*/, propagator_state_t& /*state*/,
                    const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, Args&&... /*args*/) {
    return false;
  }
};

}  // namespace Acts::detail
