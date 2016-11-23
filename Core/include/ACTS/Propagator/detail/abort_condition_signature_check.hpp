// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_ABORT_CONDITION_SIGNATURE_CHECK_HPP
#define ACTS_ABORT_CONDITION_SIGNATURE_CHECK_HPP 1

#include <type_traits>
#include "ACTS/Propagator/detail/condition_uses_result_type.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace Acts {

namespace detail {

  namespace {
    template <typename T,
              typename input,
              typename result,
              typename = decltype(std::declval<const T>().
                                  operator()(std::declval<const result&>(),
                                             std::declval<const input&>(),
                                             std::declval<double&>()))>
    std::true_type
    test_condition_with_result(int);

    template <typename, typename, typename>
    std::false_type test_condition_with_result(...);

    template <typename T,
              typename input,
              typename = decltype(std::declval<const T>().
                                  operator()(std::declval<const input&>(),
                                             std::declval<double&>()))>
    std::true_type
    test_condition_without_result(int);

    template <typename, typename>
    std::false_type test_condition_without_result(...);

    // clang-format on
    template <typename T, typename input, bool has_result = false>
    struct condition_signature_check_impl
        : decltype(test_condition_without_result<T, input>(0))
    {
    };

    template <typename T, typename input>
    struct condition_signature_check_impl<T, input, true>
        : decltype(
              test_condition_with_result<T,
                                         input,
                                         result_type_t<observer_type_t<T>>>(0))
    {
    };

    template <typename T, typename input>
    struct abort_condition_signature_check
        : condition_signature_check_impl<T,
                                         input,
                                         condition_uses_result_type<T>::value>
    {
    };
    // clang-format on
  }  // end of anonymous namespace

  template <typename T, typename input>
  constexpr bool abort_condition_signature_check_v
      = abort_condition_signature_check<T, input>::value;
}  // namespace detail

}  // namespace Acts
#endif  // ACTS_ABORT_CONDITION_SIGNATURE_CHECK_HPP
