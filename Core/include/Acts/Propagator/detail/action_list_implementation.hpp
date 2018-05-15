// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {

    /// The action caller struct, it's called with the
    /// the right result object, which is taken out
    /// from the result tuple
    template <bool has_result = true>
    struct action_caller
    {
      template <typename actor,
                typename result_t,
                typename propagator_state_t,
                typename stepper_state_t>
      static void
      action(const actor&        act,
             propagator_state_t& propState,
             stepper_state_t&    stepState,
             result_t&           r)
      {
        act(propState,
            stepState,
            r.template get<detail::result_type_t<actor>>());
      }
    };

    /// The action caller struct, without result object
    template <>
    struct action_caller<false>
    {
      template <typename actor,
                typename result_t,
                typename propagator_state_t,
                typename stepper_state_t>
      static void
      action(const actor&        act,
             propagator_state_t& propState,
             stepper_state_t&    stepState,
             result_t&)
      {
        act(propState, stepState);
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
  struct action_list_impl<first, others...>
  {
    template <typename T,
              typename result_t,
              typename propagator_state_t,
              typename stepper_state_t>
    static void
    action(const T&            obs_tuple,
           propagator_state_t& propState,
           stepper_state_t&    stepState,
           result_t&           r)
    {
      constexpr bool has_result  = has_result_type_v<first>;
      const auto&    this_action = std::get<first>(obs_tuple);
      action_caller<has_result>::action(this_action, propState, stepState, r);
      action_list_impl<others...>::action(obs_tuple, propState, stepState, r);
    }
  };

  /// The action list call implementation
  /// - it calls 'action' on the last entry of the tuple
  template <typename last>
  struct action_list_impl<last>
  {
    template <typename T,
              typename result_t,
              typename propagator_state_t,
              typename stepper_state_t>
    static void
    action(const T&            obs_tuple,
           propagator_state_t& propState,
           stepper_state_t&    stepState,
           result_t&           r)
    {
      constexpr bool has_result  = has_result_type_v<last>;
      const auto&    this_action = std::get<last>(obs_tuple);
      action_caller<has_result>::action(this_action, propState, stepState, r);
    }
  };

  /// The empty action list call implementation
  template <>
  struct action_list_impl<>
  {
    template <typename T,
              typename result_t,
              typename propagator_state_t,
              typename stepper_state_t>
    static void
    action(const T&, propagator_state_t&, stepper_state_t&, result_t&)
    {
    }
  };

}  // namespace detail
}  // namespace Acts