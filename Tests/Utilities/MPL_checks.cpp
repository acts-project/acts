// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE MPL Tests
#include <boost/mpl/equal.hpp>
#include <boost/mpl/set.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/test/included/unit_test.hpp>
#include <type_traits>
#include "ACTS/Utilities/detail/MPL/all_of.hpp"
#include "ACTS/Utilities/detail/MPL/any_of.hpp"
#include "ACTS/Utilities/detail/MPL/are_sorted.hpp"
#include "ACTS/Utilities/detail/MPL/are_within.hpp"
#include "ACTS/Utilities/detail/MPL/at_index.hpp"
#include "ACTS/Utilities/detail/MPL/boost_mpl_helper.hpp"
#include "ACTS/Utilities/detail/MPL/has_duplicates.hpp"
#include "ACTS/Utilities/detail/MPL/type_collector.hpp"

namespace bm = bm;

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_CASE(all_of_test)
  {
    using detail::all_of_v;

    static_assert(not all_of_v<true, true, false>,
                  "all_of_v<true, true, false> failed");
    static_assert(not all_of_v<false, true, true, false>,
                  "all_of_v<false, true, true, false> failed");
    static_assert(all_of_v<true, true, true>,
                  "all_of_v<true, true, true> failed");
    static_assert(all_of_v<true>, "all_of_v<true> failed");
    static_assert(not all_of_v<false>, "all_of_v<false> failed");
    static_assert(all_of_v<>, "all_of_v<> failed");
  }

  BOOST_AUTO_TEST_CASE(boost_set_merger_test)
  {
    typedef typename bm::set<float, int, char, bool>::type first;
    typedef typename bm::vector<long, int>::type               second;
    typedef typename detail::boost_set_merger_t<first, second> found;
    typedef typename bm::set<float, int, char, bool, long>::type expected;

    static_assert(std::is_same<found, expected>::value,
                  "merging sequence into bm::set failed");
  }

  template <typename... args>
  struct variadic_struct
  {
  };

  BOOST_AUTO_TEST_CASE(unpack_boost_set_as_template_test)
  {
    typedef bm::set<float, int, char>::type   boost_set;
    typedef variadic_struct<float, int, char> expected;
    typedef detail::boost_set_as_tparams_t<variadic_struct, boost_set> found;

    static_assert(std::is_same<found, expected>::value,
                  "using boost::mpl::set for variadic templates failed");
  }

  namespace {
    struct traits1
    {
      typedef int  result_type;
      typedef char observer_type;
    };

    template <bool>
    struct traits2;

    template <>
    struct traits2<false>
    {
      typedef bool  result_type;
      typedef float observer_type;
    };

    template <>
    struct traits2<true>
    {
      typedef float observer_type;
    };
  }

  BOOST_AUTO_TEST_CASE(type_collector_test)
  {
    typedef detail::type_collector_t<detail::result_type_extractor,
                                     traits1,
                                     traits2<true>,
                                     traits2<false>>
        found_results;

    typedef detail::type_collector_t<detail::observer_type_extractor,
                                     traits1,
                                     traits2<true>,
                                     traits2<false>>
        found_observers;

    typedef typename bm::set<int, bool>::type   expected_results;
    typedef typename bm::set<char, float>::type expected_observers;

    static_assert(std::is_same<found_results, expected_results>::value,
                  "collecting result types failed");
    static_assert(std::is_same<found_observers, expected_observers>::value,
                  "collecting observer types failed");
  }

  BOOST_AUTO_TEST_CASE(has_duplicates_test)
  {
    using detail::has_duplicates_v;
    static_assert(has_duplicates_v<int, float, char, int>,
                  "has_duplicates_v failed");
    static_assert(has_duplicates_v<int, int, char, float>,
                  "has_duplicates_v failed");
    static_assert(has_duplicates_v<int, char, float, float>,
                  "has_duplicates_v failed");
    static_assert(has_duplicates_v<int, char, char, float>,
                  "has_duplicates_v failed");
    static_assert(not has_duplicates_v<int, bool, char, float>,
                  "has_duplicates_v failed");
  }

  BOOST_AUTO_TEST_CASE(any_of_test)
  {
    using detail::any_of_v;

    static_assert(any_of_v<true, true, false>,
                  "any_of_v<true, true, false> failed");
    static_assert(any_of_v<false, true, true, false>,
                  "any_of_v<false, true, true, false> failed");
    static_assert(any_of_v<true, true, true>,
                  "any_of_v<true, true, true> failed");
    static_assert(not any_of_v<false, false>, "any_of_v<false, false> failed");
    static_assert(any_of_v<true>, "any_of_v<true> failed");
    static_assert(not any_of_v<false>, "any_of_v<false> failed");
    static_assert(not any_of_v<>, "any_of_v<> failed");
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
  BOOST_AUTO_TEST_CASE(are_sorted_helper_tests)
  {
    using detail::are_sorted;
    // strictly ascending
    BOOST_CHECK((are_sorted<true, true, int, -1, 3, 4, 12>::value));
    BOOST_CHECK((not are_sorted<true, true, int, -1, 13, 4>::value));
    BOOST_CHECK((not are_sorted<true, true, int, -1, 4, 4, 7>::value));
    // weakly ascending
    BOOST_CHECK((are_sorted<true, false, int, -1, 3, 4, 12>::value));
    BOOST_CHECK((not are_sorted<true, false, int, -1, 13, 4>::value));
    BOOST_CHECK((are_sorted<true, false, int, -1, 4, 4, 7>::value));
    // strictly descending
    BOOST_CHECK((are_sorted<false, true, int, 1, -3, -4, -12>::value));
    BOOST_CHECK((not are_sorted<false, true, int, 1, -13, -4>::value));
    BOOST_CHECK((not are_sorted<false, true, int, 1, -4, -4>::value));
    // weakly descending
    BOOST_CHECK((are_sorted<false, false, int, 1, -3, -4, -12>::value));
    BOOST_CHECK((not are_sorted<false, false, int, -1, -13, -4>::value));
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
  BOOST_AUTO_TEST_CASE(are_within_helper_tests)
  {
    using detail::are_within;
    BOOST_CHECK((are_within<int, 0, 10, 1, 3, 7, 2>::value));
    BOOST_CHECK((are_within<int, 0, 10, 1, 3, 0, 2>::value));
    BOOST_CHECK((not are_within<int, 0, 10, -1, 3, 7, 2>::value));
    BOOST_CHECK((not are_within<int, 0, 10, -1, 3, 7, -2>::value));
    BOOST_CHECK((not are_within<int, 0, 10, 1, 3, 17, 2>::value));
    BOOST_CHECK((not are_within<int, 0, 10, 1, 3, 17, 12>::value));
    BOOST_CHECK((not are_within<int, 0, 10, 1, 10>::value));
    BOOST_CHECK((not are_within<int, 0, 10, 1, -2, 10, 14>::value));
  }

  /**
   * @brief Unit test for Acts::details::at_index helper
   */
  BOOST_AUTO_TEST_CASE(at_index_helper_tests)
  {
    using detail::at_index;
    BOOST_CHECK((at_index<int, 0, 10, 1, 3, 7, 2>::value == 10));
    BOOST_CHECK((at_index<int, 1, 10, 1, 3, 7, 2>::value == 1));
    BOOST_CHECK((at_index<int, 2, 10, 1, 3, 7, 2>::value == 3));
    BOOST_CHECK((at_index<int, 3, 10, 1, 3, 7, 2>::value == 7));
    BOOST_CHECK((at_index<int, 4, 10, 1, 3, 7, 2>::value == 2));
  }
}  // namespace Test

}  // namespace Acts
