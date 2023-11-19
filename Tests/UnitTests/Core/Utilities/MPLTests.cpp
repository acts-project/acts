// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/any_of.hpp"
#include "Acts/Utilities/detail/MPL/are_sorted.hpp"
#include "Acts/Utilities/detail/MPL/are_within.hpp"
#include "Acts/Utilities/detail/MPL/at_index.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

#include <tuple>
#include <type_traits>

#include <boost/hana.hpp>
#include <boost/hana/core/to.hpp>
#include <boost/hana/equal.hpp>
#include <boost/hana/ext/std/tuple.hpp>
#include <boost/hana/integral_constant.hpp>
#include <boost/hana/set.hpp>
#include <boost/hana/transform.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/type.hpp>
#include <boost/hana/union.hpp>

namespace hana = boost::hana;
namespace Acts {
namespace detail {
template <bool ascending, bool strict, typename T, T... values>
struct are_sorted;
template <typename T, T MIN, T MAX, T... values>
struct are_within;
template <typename T, std::size_t index, T... values>
struct at_index;
}  // namespace detail

namespace Test {

BOOST_AUTO_TEST_CASE(all_of_test) {
  using detail::all_of_v;

  static_assert(!all_of_v<true, true, false>,
                "all_of_v<true, true, false> failed");
  static_assert(!all_of_v<false, true, true, false>,
                "all_of_v<false, true, true, false> failed");
  static_assert(all_of_v<true, true, true>,
                "all_of_v<true, true, true> failed");
  static_assert(all_of_v<true>, "all_of_v<true> failed");
  static_assert(!all_of_v<false>, "all_of_v<false> failed");
  static_assert(all_of_v<>, "all_of_v<> failed");
}

BOOST_AUTO_TEST_CASE(hana_set_union_test) {
  // using first    = typename bm::set<float, int, char, bool>::type;
  constexpr auto first = hana::make_set(hana::type_c<float>, hana::type_c<int>,
                                        hana::type_c<char>, hana::type_c<bool>);
  // using second   = typename bm::vector<long, int>::type;
  constexpr auto second = hana::make_set(hana::type_c<long>, hana::type_c<int>);
  constexpr auto found = hana::union_(first, second);
  // using found    = typename detail::boost_set_merger_t<first, second>;
  // using expected = typename bm::set<float, int, char, bool, long>::type;
  constexpr auto expected =
      hana::make_set(hana::type_c<float>, hana::type_c<int>, hana::type_c<char>,
                     hana::type_c<bool>, hana::type_c<long>);

  static_assert(found == expected, "union of hana::sets failed");
}

BOOST_AUTO_TEST_CASE(hana_set_to_tuple_test) {
  constexpr auto a_set = hana::make_set(hana::type_c<float>, hana::type_c<int>,
                                        hana::type_c<char>, hana::type_c<bool>);
  constexpr auto h_tuple =
      hana::make_tuple(hana::type_c<float>, hana::type_c<int>,
                       hana::type_c<char>, hana::type_c<bool>);

  static_assert(hana::to<hana::tuple_tag>(a_set) == h_tuple, "not equal");

  // using std_tuple = decltype(hana::unpack(a_set,
  // hana::template_<std::tuple>))::type; using expected = std::tuple<float,
  // int, char>; static_assert(std::is_same<std_tuple, expected>::value,
  // "using
  // boost::mpl::set for variadic templates failed");
}

template <typename... args>
struct variadic_struct {
  using tuple = std::tuple<args...>;
};

BOOST_AUTO_TEST_CASE(unpack_boost_set_as_template_test) {
  constexpr auto hana_set = hana::make_set(
      hana::type_c<float>, hana::type_c<int>, hana::type_c<char>);
  using found =
      decltype(hana::unpack(hana_set, hana::template_<variadic_struct>))::type;

  using expected = variadic_struct<float, int, char>;

  static_assert(std::is_same<found, expected>::value,
                "using boost::mpl::set for variadic templates failed");

  static_assert(
      std::is_same<expected::tuple, std::tuple<float, int, char>>::value,
      "not equal");
}

namespace {
struct traits1 {
  using result_type = int;
  using action_type = char;
};

template <bool>
struct traits2;

template <>
struct traits2<false> {
  using result_type = bool;
  using action_type = float;
};

template <>
struct traits2<true> {
  using action_type = float;
};
}  // namespace

template <typename... Args>
struct tuple_helper {
  using tuple = std::tuple<Args...>;
};

BOOST_AUTO_TEST_CASE(type_collector_test) {
  // test some predicates
  static_assert(detail::has_result_type_v<traits1>, "Did not find result type");
  static_assert(detail::has_result_type_v<traits2<false>>,
                "Did not find result type");
  static_assert(!detail::has_result_type_v<traits2<true>>,
                "Did find result type");

  static_assert(detail::has_action_type_v<traits1>, "Did not find action type");
  static_assert(detail::has_action_type_v<traits2<false>>,
                "Did not find action type");
  static_assert(detail::has_action_type_v<traits2<true>>,
                "Did not find action type");

  constexpr auto found_results =
      detail::type_collector_t<detail::result_type_extractor, traits1,
                               traits2<true>, traits2<false>>;
  constexpr auto expected_results =
      hana::make_set(hana::type_c<int>, hana::type_c<bool>);
  static_assert(found_results == expected_results,
                "Didn't find expected results");

  // check unpack
  using found_results_tuple = decltype(hana::unpack(
      found_results, hana::template_<tuple_helper>))::type::tuple;
  using expected_results_tuple = std::tuple<int, bool>;
  static_assert(
      std::is_same<found_results_tuple, expected_results_tuple>::value,
      "Unpacked results tuple not correct");

  constexpr auto found_actions =
      detail::type_collector_t<detail::action_type_extractor, traits1,
                               traits2<true>, traits2<false>>;
  constexpr auto expected_actions =
      hana::make_set(hana::type_c<char>, hana::type_c<float>);
  static_assert(found_actions == expected_actions,
                "Didn't find expected actions");

  // check unpack
  using found_actions_tuple = decltype(hana::unpack(
      found_actions, hana::template_<tuple_helper>))::type::tuple;
  using expected_actions_tuple = std::tuple<char, float>;
  static_assert(
      std::is_same<found_actions_tuple, expected_actions_tuple>::value,
      "Unpacked actions tuple not correct");
}

BOOST_AUTO_TEST_CASE(has_duplicates_test) {
  using detail::has_duplicates_v;
  static_assert(has_duplicates_v<int, float, char, int>,
                "has_duplicates_v failed");
  static_assert(has_duplicates_v<int, int, char, float>,
                "has_duplicates_v failed");
  static_assert(has_duplicates_v<int, char, float, float>,
                "has_duplicates_v failed");
  static_assert(has_duplicates_v<int, char, char, float>,
                "has_duplicates_v failed");
  static_assert(!has_duplicates_v<int, bool, char, float>,
                "has_duplicates_v failed");
}

BOOST_AUTO_TEST_CASE(any_of_test) {
  using detail::any_of_v;

  static_assert(any_of_v<true, true, false>,
                "any_of_v<true, true, false> failed");
  static_assert(any_of_v<false, true, true, false>,
                "any_of_v<false, true, true, false> failed");
  static_assert(any_of_v<true, true, true>,
                "any_of_v<true, true, true> failed");
  static_assert(!any_of_v<false, false>, "any_of_v<false, false> failed");
  static_assert(any_of_v<true>, "any_of_v<true> failed");
  static_assert(!any_of_v<false>, "any_of_v<false> failed");
  static_assert(!any_of_v<>, "any_of_v<> failed");
}

/**
 * @brief Unit test for Acts::anonymous_namespace{ParameterSet.h}::are_sorted
 * helper
 *
 * The test checks for correct behavior in the following cases (always using
 * @c int
 * as value type):
 * -# test: ordered strictly ascending, input: ordered strictly ascending
 * -# test: ordered strictly ascending, input: unordered
 * -# test: ordered strictly ascending, input: ordered weakly ascending
 * -# test: ordered weakly ascending, input: ordered strictly ascending
 * -# test: ordered weakly ascending, input: unordered
 * -# test: ordered weakly ascending, input: ordered weakly ascending
 * -# test: ordered strictly descending, input: ordered strictly descending
 * -# test: ordered strictly descending, input: unordered
 * -# test: ordered strictly descending, input: ordered weakly descending
 * -# test: ordered weakly descending, input: ordered strictly descending
 * -# test: ordered weakly descending, input: unordered
 * -# test: ordered weakly descending, input: ordered weakly descending
 */
BOOST_AUTO_TEST_CASE(are_sorted_helper_tests) {
  using detail::are_sorted;
  // strictly ascending
  BOOST_CHECK((are_sorted<true, true, int, -1, 3, 4, 12>::value));
  BOOST_CHECK((!are_sorted<true, true, int, -1, 13, 4>::value));
  BOOST_CHECK((!are_sorted<true, true, int, -1, 4, 4, 7>::value));
  // weakly ascending
  BOOST_CHECK((are_sorted<true, false, int, -1, 3, 4, 12>::value));
  BOOST_CHECK((!are_sorted<true, false, int, -1, 13, 4>::value));
  BOOST_CHECK((are_sorted<true, false, int, -1, 4, 4, 7>::value));
  // strictly descending
  BOOST_CHECK((are_sorted<false, true, int, 1, -3, -4, -12>::value));
  BOOST_CHECK((!are_sorted<false, true, int, 1, -13, -4>::value));
  BOOST_CHECK((!are_sorted<false, true, int, 1, -4, -4>::value));
  // weakly descending
  BOOST_CHECK((are_sorted<false, false, int, 1, -3, -4, -12>::value));
  BOOST_CHECK((!are_sorted<false, false, int, -1, -13, -4>::value));
  BOOST_CHECK((are_sorted<false, false, int, -1, -4, -4, -7>::value));
}

/**
 * @brief Unit test for Acts::anonymous_namespace{ParameterSet.h}::are_within
 * helper
 *
 * The test checks for correct behavior in the following cases (always using
 * @c int
 * as value type):
 * -# all values within (MIN,MAX)
 * -# all values within [MIN,MAX)
 * -# one value < MIN
 * -# multiple values < MIN
 * -# one value > MAX
 * -# multiple values > Max
 * -# one value == MAX
 * -# contains values < MIN and >= MAX
 */
BOOST_AUTO_TEST_CASE(are_within_helper_tests) {
  using detail::are_within;
  BOOST_CHECK((are_within<int, 0, 10, 1, 3, 7, 2>::value));
  BOOST_CHECK((are_within<int, 0, 10, 1, 3, 0, 2>::value));
  BOOST_CHECK((!are_within<int, 0, 10, -1, 3, 7, 2>::value));
  BOOST_CHECK((!are_within<int, 0, 10, -1, 3, 7, -2>::value));
  BOOST_CHECK((!are_within<int, 0, 10, 1, 3, 17, 2>::value));
  BOOST_CHECK((!are_within<int, 0, 10, 1, 3, 17, 12>::value));
  BOOST_CHECK((!are_within<int, 0, 10, 1, 10>::value));
  BOOST_CHECK((!are_within<int, 0, 10, 1, -2, 10, 14>::value));
}

/**
 * @brief Unit test for Acts::details::at_index helper
 */
BOOST_AUTO_TEST_CASE(at_index_helper_tests) {
  using detail::at_index;
  BOOST_CHECK_EQUAL((at_index<int, 0, 10, 1, 3, 7, 2>::value), 10);
  BOOST_CHECK_EQUAL((at_index<int, 1, 10, 1, 3, 7, 2>::value), 1);
  BOOST_CHECK_EQUAL((at_index<int, 2, 10, 1, 3, 7, 2>::value), 3);
  BOOST_CHECK_EQUAL((at_index<int, 3, 10, 1, 3, 7, 2>::value), 7);
  BOOST_CHECK_EQUAL((at_index<int, 4, 10, 1, 3, 7, 2>::value), 2);
}
}  // namespace Test

}  // namespace Acts
