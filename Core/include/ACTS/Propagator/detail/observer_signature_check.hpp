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
              typename input,
              typename result,
              typename = decltype(std::declval<T>().
                                  operator()(std::declval<const input&>(),
                                             std::declval<result&>()))>
    std::true_type
    test_observer_with_result(int);

    template <typename, typename, typename>
    std::false_type
    test_observer_with_result(...);

    template <typename T,
              typename input,
              typename = decltype(std::declval<T>().
                                  operator()(std::declval<const input&>()))>
    std::true_type
    test_observer_without_result(int);

    template <typename, typename>
    std::false_type
    test_observer_without_result(...);

    template <typename T, typename input, bool has_result = false>
    struct observer_signature_check_impl
        : decltype(test_observer_without_result<T, input>(0))
    {
    };

    template <typename T, typename input>
    struct observer_signature_check_impl<T, input, true>
        : decltype(
              test_observer_with_result<T, input, detail::result_type_t<T>>(0))
    {
    };

    template <typename T, typename input>
    struct observer_signature_check
        : observer_signature_check_impl<T, input, detail::has_result_type_v<T>>
    {
    };
  }  // end of anonymous namespace

  template <typename T, typename input>
  constexpr bool observer_signature_check_v
      = observer_signature_check<T, input>::value;
}  // namespace detail

}  // namespace Acts
#endif  // ACTS_OBSERVER_SIGNATURE_CHECK_HPP
