// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_ABORT_LIST_IMPLEMENTATION_HPP
#define ACTS_ABORT_LIST_IMPLEMENTATION_HPP 1

#include <algorithm>
#include "ACTS/Propagator/detail/condition_uses_result_type.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {
    template <bool has_result = true>
    struct condition_caller
    {
      template <typename condition, typename result, typename input>
      static bool
      check(const condition& c,
            const result&    r,
            input&           cache,
            double&          stepMax)
      {
        typedef observer_type_t<condition>   observer_type;
        typedef result_type_t<observer_type> result_type;

        return c(r.template get<result_type>(), cache, stepMax);
      }
    };

    template <>
    struct condition_caller<false>
    {
      template <typename condition, typename result, typename input>
      static bool
      check(const condition& c, const result&, input& cache, double& stepMax)
      {
        return c(cache, stepMax);
      }
    };
  }  // end of anonymous namespace

  template <typename... conditions>
  struct abort_list_impl;

  template <typename first, typename... others>
  struct abort_list_impl<first, others...>
  {
    template <typename T, typename result, typename input>
    static bool
    check(const T&      conditions_tuple,
          const result& r,
          input&        cache,
          double&       stepMax)
    {
      // get the right helper for calling the abort condition
      constexpr bool has_result = condition_uses_result_type<first>::value;
      typedef condition_caller<has_result> caller_type;

      // get the cache abort condition
      const auto& this_condition = std::get<first>(conditions_tuple);

      // maximum step sizes from remaining abort conditions
      double other_stepMax = stepMax;

      // - check abort conditions recursively
      // - make use of short-circuit evaluation
      // -> skip remaining conditions if this abort condition evaluates to true
      bool abort = caller_type::check(this_condition, r, cache, stepMax)
          || abort_list_impl<others...>::check(
                       conditions_tuple, r, cache, other_stepMax);

      // set remaining step size
      stepMax = std::min(stepMax, other_stepMax);

      return abort;
    }
  };

  template <typename last>
  struct abort_list_impl<last>
  {
    template <typename T, typename result, typename input>
    static bool
    check(const T&      conditions_tuple,
          const result& r,
          input&        cache,
          double&       stepMax)
    {
      // get the right helper for calling the abort condition
      constexpr bool has_result     = condition_uses_result_type<last>::value;
      const auto&    this_condition = std::get<last>(conditions_tuple);

      return condition_caller<has_result>::check(
          this_condition, r, cache, stepMax);
    }
  };

  template <>
  struct abort_list_impl<>
  {
    template <typename T, typename result, typename input>
    static bool
    check(const T&, const result&, input&, double&)
    {
      return false;
    }
  };

}  // namespace

}  // namespace Acts
#endif  // ACTS_ABORT_LIST_IMPLEMENTATION_HPP
