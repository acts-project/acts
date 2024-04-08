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
 * @brief check whether integral values are sorted
 *
 * @tparam ascending boolean flag to check for ascending order (@c true) or
 * descending order (@c false)
 * @tparam strict boolean flag whether strict ordering is required
 * @tparam T integral type of values whose order should be checked
 * @tparam values template parameter pack containing the list of values
 *
 * @test Unit tests are implemented \link
 * Acts::Test::BOOST_AUTO_TEST_CASE(are_sorted_helper_tests) here\endlink.
 *
 * @return `are_sorted<asc,strict,T,values...>::value` is @c true if the given
 * values are properly sorted,
 *         otherwise @c false
 */
template <bool ascending, bool strict, typename T, T... values>
struct are_sorted;

/// @cond
// one value is always sorted
template <bool ascending, bool strict, typename T, T v>
struct are_sorted<ascending, strict, T, v> {
  enum { value = true };
};

// strict, ascending ordering
template <typename T, T a, T b, T... N>
struct are_sorted<true, true, T, a, b, N...> {
  enum { value = ((a < b) && are_sorted<true, true, T, b, N...>::value) };
};

// weak, ascending ordering
template <typename T, T a, T b, T... N>
struct are_sorted<true, false, T, a, b, N...> {
  enum { value = (a <= b && are_sorted<true, false, T, b, N...>::value) };
};

// strict, descending ordering
template <typename T, T a, T b, T... N>
struct are_sorted<false, true, T, a, b, N...> {
  enum { value = (a > b && are_sorted<false, true, T, b, N...>::value) };
};

// weak, descending ordering
template <typename T, T a, T b, T... N>
struct are_sorted<false, false, T, a, b, N...> {
  enum { value = (a >= b && are_sorted<false, false, T, b, N...>::value) };
};
/// @endcond
}  // namespace Acts::detail
