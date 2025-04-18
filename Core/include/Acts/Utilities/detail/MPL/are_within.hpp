// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
namespace Acts::detail {
/**
 * @brief check whether integral values are within a given range
 *
 * @tparam T integral type of values whose range should be checked
 * @tparam MIN lower accepted bound of values (inclusive)
 * @tparam MAX upper accepted bound of values (exclusive)
 * @tparam values template parameter pack containing the list of values
 *
 * @test Unit tests are implemented \link
 * Acts::Test::BOOST_AUTO_TEST_CASE(are_within_helper_tests) here\endlink.
 *
 * @return `are_within<T,MIN,MAX,values...>::value` is @c true if all given
 * values are within the
 *          interval [MIN,MAX), otherwise @c false
 */
template <typename T, T MIN, T MAX, T... values>
struct are_within;

/// @cond
// check last parameter
template <typename T, T MIN, T MAX, T a>
struct are_within<T, MIN, MAX, a> {
  enum { value = (a >= MIN && a < MAX) };
};

// recursive check
template <typename T, T MIN, T MAX, T a, T... others>
struct are_within<T, MIN, MAX, a, others...> {
  enum {
    value =
        ((a >= MIN) && (a < MAX) && are_within<T, MIN, MAX, others...>::value)
  };
};
/// @endcond
}  // namespace Acts::detail
