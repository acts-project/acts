// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_OBSERVER_SIGNATURE_CHECK_HPP
#define ACTS_OBSERVER_SIGNATURE_CHECK_HPP 1

#include <type_traits>
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

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
    test_observer_with_result(int);

    template <typename, typename, typename>
    std::false_type
    test_observer_with_result(...);

    template <typename T,
              typename cache_t,
              typename = decltype(std::declval<T>().
                                  operator()(std::declval<cache_t&>()))>
    std::true_type
    test_observer_without_result(int);

    template <typename, typename>
    std::false_type
    test_observer_without_result(...);

    template <typename T, typename cache_t, bool has_result = false>
    struct observer_signature_check_impl
        : decltype(test_observer_without_result<T, cache_t>(0))
    {
    };

    template <typename T, typename cache_t>
    struct observer_signature_check_impl<T, cache_t, true>
        : decltype(
              test_observer_with_result<T, cache_t, detail::result_type_t<T>>(
                  0))
    {
    };

    template <typename T, typename cache_t>
    struct observer_signature_check
        : observer_signature_check_impl<T,
                                        cache_t,
                                        detail::has_result_type_v<T>>
    {
    };
  }  // end of anonymous namespace

  template <typename T, typename cache_t>
  constexpr bool observer_signature_check_v
      = observer_signature_check<T, cache_t>::value;
}  // namespace detail

}  // namespace Acts
#endif  // ACTS_OBSERVER_SIGNATURE_CHECK_HPP
