// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/JoinStrings.hpp"

#include <array>
#include <list>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(EmptyRange) {
  std::vector<std::string> empty;
  BOOST_CHECK_EQUAL(joinStrings(empty, ","), "");

  std::vector<int> emptyInts;
  BOOST_CHECK_EQUAL(joinStrings(emptyInts, ","), "");

  // Test with custom format
  BOOST_CHECK_EQUAL(joinStrings(emptyInts, ",", "{:02d}"), "");
}

BOOST_AUTO_TEST_CASE(SingleElement) {
  std::vector<std::string> single = {"hello"};
  BOOST_CHECK_EQUAL(joinStrings(single, ","), "hello");

  std::vector<int> singleInt = {42};
  BOOST_CHECK_EQUAL(joinStrings(singleInt, ","), "42");

  // Test with custom format
  BOOST_CHECK_EQUAL(joinStrings(singleInt, ",", "{:04d}"), "0042");
}

BOOST_AUTO_TEST_CASE(MultipleElements) {
  std::vector<std::string> strings = {"apple", "banana", "cherry"};
  BOOST_CHECK_EQUAL(joinStrings(strings, ","), "apple,banana,cherry");
  BOOST_CHECK_EQUAL(joinStrings(strings, ", "), "apple, banana, cherry");
  BOOST_CHECK_EQUAL(joinStrings(strings, " | "), "apple | banana | cherry");
}

BOOST_AUTO_TEST_CASE(IntegerTypes) {
  std::vector<int> ints = {1, 2, 3, 4, 5};
  BOOST_CHECK_EQUAL(joinStrings(ints, ","), "1,2,3,4,5");
  BOOST_CHECK_EQUAL(joinStrings(ints, " - "), "1 - 2 - 3 - 4 - 5");

  // Test with custom format
  BOOST_CHECK_EQUAL(joinStrings(ints, ",", "{:02d}"), "01,02,03,04,05");
  BOOST_CHECK_EQUAL(joinStrings(ints, ",", "{:+d}"), "+1,+2,+3,+4,+5");
}

BOOST_AUTO_TEST_CASE(FloatingPointTypes) {
  std::vector<double> doubles = {1.5, 2.7, 3.14};
  BOOST_CHECK_EQUAL(joinStrings(doubles, ","), "1.5,2.7,3.14");

  // Test with custom format
  BOOST_CHECK_EQUAL(joinStrings(doubles, ",", "{:.1f}"), "1.5,2.7,3.1");
  BOOST_CHECK_EQUAL(joinStrings(doubles, ",", "{:.2f}"), "1.50,2.70,3.14");
}

BOOST_AUTO_TEST_CASE(MixedNumericTypes) {
  std::vector<long> longs = {1000L, 2000L, 3000L};
  BOOST_CHECK_EQUAL(joinStrings(longs, " | "), "1000 | 2000 | 3000");

  std::vector<unsigned int> uints = {10u, 20u, 30u};
  BOOST_CHECK_EQUAL(joinStrings(uints, " -> "), "10 -> 20 -> 30");
}

BOOST_AUTO_TEST_CASE(DifferentContainerTypes) {
  // std::array
  std::array<std::string, 3> arr = {"x", "y", "z"};
  BOOST_CHECK_EQUAL(joinStrings(arr, "-"), "x-y-z");

  // std::list
  std::list<int> lst = {7, 8, 9};
  BOOST_CHECK_EQUAL(joinStrings(lst, " "), "7 8 9");

  // C-style array via span-like interface
  int cArray[] = {10, 11, 12};
  std::vector<int> vecFromArray(std::begin(cArray), std::end(cArray));
  BOOST_CHECK_EQUAL(joinStrings(vecFromArray, "/"), "10/11/12");
}

BOOST_AUTO_TEST_CASE(StringViewCompatible) {
  // Test with string literals and const char*
  std::vector<const char*> cstrings = {"hello", "world", "test"};
  BOOST_CHECK_EQUAL(joinStrings(cstrings, " "), "hello world test");

  // Test with std::string_view
  std::vector<std::string_view> stringViews = {"foo", "bar", "baz"};
  BOOST_CHECK_EQUAL(joinStrings(stringViews, "::"), "foo::bar::baz");
}

BOOST_AUTO_TEST_CASE(SpecialDelimiters) {
  std::vector<std::string> strings = {"a", "b", "c"};

  // Empty delimiter
  BOOST_CHECK_EQUAL(joinStrings(strings, ""), "abc");

  // Multi-character delimiters
  BOOST_CHECK_EQUAL(joinStrings(strings, " and "), "a and b and c");
  BOOST_CHECK_EQUAL(joinStrings(strings, " -> "), "a -> b -> c");

  // Special characters
  BOOST_CHECK_EQUAL(joinStrings(strings, "\n"), "a\nb\nc");
  BOOST_CHECK_EQUAL(joinStrings(strings, "\t"), "a\tb\tc");
}

BOOST_AUTO_TEST_CASE(EmptyStrings) {
  std::vector<std::string> withEmpty = {"", "middle", ""};
  BOOST_CHECK_EQUAL(joinStrings(withEmpty, ","), ",middle,");

  std::vector<std::string> allEmpty = {"", "", ""};
  BOOST_CHECK_EQUAL(joinStrings(allEmpty, "-"), "--");
}

BOOST_AUTO_TEST_CASE(CustomFormatStrings) {
  std::vector<int> numbers = {1, 2, 3};

  // Hexadecimal
  BOOST_CHECK_EQUAL(joinStrings(numbers, ",", "{:x}"), "1,2,3");
  BOOST_CHECK_EQUAL(joinStrings(numbers, ",", "{:X}"), "1,2,3");

  // With prefix
  BOOST_CHECK_EQUAL(joinStrings(numbers, ",", "#{:02d}"), "#01,#02,#03");

  std::vector<double> values = {1.234, 5.678};
  BOOST_CHECK_EQUAL(joinStrings(values, " | ", "{:.1e}"), "1.2e+00 | 5.7e+00");
}

BOOST_AUTO_TEST_CASE(LargerRanges) {
  // Test with larger ranges to ensure performance is reasonable
  std::vector<int> large(100);
  std::iota(large.begin(), large.end(), 1);

  std::string result = joinStrings(large, ",");

  // Check start and end of result
  BOOST_CHECK(result.starts_with("1,2,3"));
  BOOST_CHECK(result.ends_with("98,99,100"));

  // Check total length is reasonable (should contain all numbers and
  // delimiters)
  BOOST_CHECK(result.length() > 100);  // At least 100 digits + 99 commas
}

BOOST_AUTO_TEST_CASE(BooleanTypes) {
  // Use std::array instead of std::vector<bool> due to proxy reference issues
  std::array<bool, 3> bools = {true, false, true};
  BOOST_CHECK_EQUAL(joinStrings(bools, ","), "true,false,true");

  // With custom format to get numeric representation
  BOOST_CHECK_EQUAL(joinStrings(bools, " | ", "{:d}"), "1 | 0 | 1");
}

BOOST_AUTO_TEST_CASE(CharTypes) {
  std::vector<char> chars = {'a', 'b', 'c'};
  BOOST_CHECK_EQUAL(joinStrings(chars, "-"), "a-b-c");

  // With custom format for ASCII values
  BOOST_CHECK_EQUAL(joinStrings(chars, ",", "{:d}"), "97,98,99");
}

BOOST_AUTO_TEST_CASE(ViewsTransform) {
  // Test with std::views::transform - similar to real usage in GraphViz.cpp
  std::vector<int> numbers = {1, 2, 3, 4};

  auto squared = std::views::transform(numbers, [](int n) { return n * n; });
  BOOST_CHECK_EQUAL(joinStrings(squared, ","), "1,4,9,16");

  // Test with string transformation
  std::vector<std::string> words = {"hello", "world", "test"};
  auto uppercased = std::views::transform(words, [](const std::string& s) {
    std::string result;
    std::transform(s.begin(), s.end(), std::back_inserter(result),
                   [](char c) { return std::toupper(c); });
    return result;
  });
  BOOST_CHECK_EQUAL(joinStrings(uppercased, " "), "HELLO WORLD TEST");

  // Test with custom format on transformed view
  auto prefixed = std::views::transform(numbers, [](int n) { return n * 10; });
  BOOST_CHECK_EQUAL(joinStrings(prefixed, " | ", "#{:02d}"),
                    "#10 | #20 | #30 | #40");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
