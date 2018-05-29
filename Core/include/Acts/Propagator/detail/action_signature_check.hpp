// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <type_traits>
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {
    template <typename T,
              typename cache_t,
              typename result_t,
              typename = decltype(std::declval<T>().
                                  operator()(std::declval<cache_t&>(),
                                             std::declval<result_t&>()))>
    std::true_type
    test_action_with_result(int);

    template <typename, typename, typename>
    std::false_type
    test_action_with_result(...);

    template <typename T,
              typename cache_t,
              typename = decltype(std::declval<T>().
                                  operator()(std::declval<cache_t&>()))>
    std::true_type
    test_action_without_result(int);

    template <typename, typename>
    std::false_type
    test_action_without_result(...);

    template <typename T, typename cache_t, bool has_result = false>
    struct action_signature_check_impl
        : decltype(test_action_without_result<T, cache_t>(0))
    {
    };

    template <typename T, typename cache_t>
    struct action_signature_check_impl<T, cache_t, true>
        : decltype(
              test_action_with_result<T, cache_t, detail::result_type_t<T>>(0))
    {
    };

    template <typename T, typename cache_t>
    struct action_signature_check
        : action_signature_check_impl<T, cache_t, detail::has_result_type_v<T>>
    {
    };
  }  // end of anonymous namespace

  template <typename T, typename cache_t>
  constexpr bool action_signature_check_v
      = action_signature_check<T, cache_t>::value;
}  // namespace detail

}  // namespace Acts