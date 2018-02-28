// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_OBSERVER_LIST_IMPLEMENTATION_HPP
#define ACTS_OBSERVER_LIST_IMPLEMENTATION_HPP

#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {

    template <bool has_result = true>
    struct observer_caller
    {
      template <typename observer, typename result, typename input>
      static void
      observe(const observer& obs, const input& cache, result& r)
      {
        obs(cache, r.template get<detail::result_type_t<observer>>());
      }
    };

    template <>
    struct observer_caller<false>
    {
      template <typename observer, typename result, typename input>
      static void
      observe(const observer& obs, const input& cache, result&)
      {
        obs(cache);
      }
    };
  }  // end of anonymous namespace

  template <typename... observers>
  struct observer_list_impl;

  template <typename first, typename... others>
  struct observer_list_impl<first, others...>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T& obs_tuple, const input& cache, result& r)
    {
      constexpr bool has_result    = has_result_type_v<first>;
      const auto&    this_observer = std::get<first>(obs_tuple);
      observer_caller<has_result>::observe(this_observer, cache, r);
      observer_list_impl<others...>::observe(obs_tuple, cache, r);
    }
  };

  template <typename last>
  struct observer_list_impl<last>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T& obs_tuple, const input& cache, result& r)
    {
      constexpr bool has_result    = has_result_type_v<last>;
      const auto&    this_observer = std::get<last>(obs_tuple);
      observer_caller<has_result>::observe(this_observer, cache, r);
    }
  };

  template <>
  struct observer_list_impl<>
  {
    template <typename T, typename result, typename input>
    static void
    observe(const T&, const input&, result&)
    {
    }
  };
}  // namespace detail

}  // namespace Acts
#endif  // ACTS_OBSERVER_LIST_IMPLEMENTATION_HPP
