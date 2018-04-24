// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_ABORT_LIST_IMPLEMENTATION_HPP
#define ACTS_ABORT_LIST_IMPLEMENTATION_HPP

#include <algorithm>
#include "Acts/Propagator/detail/condition_uses_result_type.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {
    /// This is the caller used if the condition uses the result
    /// - and that result type exists
    template <bool has_result = true>
    struct condition_caller
    {
      template <typename condition, typename result_t, typename cache_t>
      static bool
      check(const condition& c, const result_t& r, cache_t& cache)
      {
        typedef action_type_t<condition>   action_type;
        typedef result_type_t<action_type> result_type;

        return c(r.template get<result_type>(), cache);
      }
    };

    /// This is the caller used if the condition only uses the cache
    /// it has no access to the result
    template <>
    struct condition_caller<false>
    {
      template <typename condition, typename result_t, typename cache_t>
      static bool
      check(const condition& c, const result_t&, cache_t& cache)
      {
        return c(cache);
      }
    };
  }  // end of anonymous namespace

  template <typename... conditions>
  struct abort_list_impl;

  template <typename first, typename... others>
  struct abort_list_impl<first, others...>
  {
    template <typename T, typename result_t, typename cache_t>
    static bool
    check(const T& conditions_tuple, const result_t& r, cache_t& cache)
    {

      // get the right helper for calling the abort condition
      constexpr bool has_result = condition_uses_result_type<first>::value;
      typedef condition_caller<has_result> caller_type;

      // get the cache abort condition
      const auto& this_condition = std::get<first>(conditions_tuple);

      // - check abort conditions recursively
      // - make use of short-circuit evaluation
      // -> skip remaining conditions if this abort condition evaluates to true
      bool abort = caller_type::check(this_condition, r, cache)
          || abort_list_impl<others...>::check(conditions_tuple, r, cache);

      return abort;
    }
  };

  template <typename last>
  struct abort_list_impl<last>
  {
    template <typename T, typename result_t, typename cache_t>
    static bool
    check(const T& conditions_tuple, const result_t& r, cache_t& cache)
    {
      // get the right helper for calling the abort condition
      constexpr bool has_result     = condition_uses_result_type<last>::value;
      const auto&    this_condition = std::get<last>(conditions_tuple);

      return condition_caller<has_result>::check(this_condition, r, cache);
    }
  };

  template <>
  struct abort_list_impl<>
  {
    template <typename T, typename result_t, typename cache_t>
    static bool
    check(const T&, const result_t&, cache_t&)
    {
      return false;
    }
  };

}  // namespace

}  // namespace Acts
#endif  // ACTS_ABORT_LIST_IMPLEMENTATION_HPP
