// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_ACTION_LIST_IMPLEMENTATION_HPP
#define ACTS_ACTION_LIST_IMPLEMENTATION_HPP

#include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {

    template <bool has_result = true>
    struct action_caller
    {
      template <typename actor, typename result_t, typename cache_t>
      static void
      action(const actor& act, cache_t& cache, result_t& r)
      {
        act(cache, r.template get<detail::result_type_t<actor>>());
      }
    };

    template <>
    struct action_caller<false>
    {
      template <typename actor, typename result_t, typename cache_t>
      static void
      action(const actor& act, cache_t& cache, result_t&)
      {
        act(cache);
      }
    };
  }  // end of anonymous namespace

  template <typename... actors>
  struct action_list_impl;

  template <typename first, typename... others>
  struct action_list_impl<first, others...>
  {
    template <typename T, typename result_t, typename cache_t>
    static void
    action(const T& obs_tuple, cache_t& cache, result_t& r)
    {
      constexpr bool has_result  = has_result_type_v<first>;
      const auto&    this_action = std::get<first>(obs_tuple);
      action_caller<has_result>::action(this_action, cache, r);
      action_list_impl<others...>::action(obs_tuple, cache, r);
    }
  };

  template <typename last>
  struct action_list_impl<last>
  {
    template <typename T, typename result_t, typename cache_t>
    static void
    action(const T& obs_tuple, cache_t& cache, result_t& r)
    {
      constexpr bool has_result  = has_result_type_v<last>;
      const auto&    this_action = std::get<last>(obs_tuple);
      action_caller<has_result>::action(this_action, cache, r);
    }
  };

  template <>
  struct action_list_impl<>
  {
    template <typename T, typename result_t, typename cache_t>
    static void
    action(const T&, cache_t&, result_t&)
    {
    }
  };

}  // namespace detail
}  // namespace Acts

#endif  // ACTS_ACTION_LIST_IMPLEMENTATION_HPP
