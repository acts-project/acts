// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for dfe::FlatMap

#include <boost/test/unit_test.hpp>

#include <cstdint>
#include <string>

#include "dfe/dfe_flat.hpp"

BOOST_AUTO_TEST_CASE(flatmap_int_int) {
  dfe::FlatMap<int, int> m;

  m.emplace(2, 12);
  m.emplace(-2, 23);
  m.emplace(54321, -2);
  m.emplace(INT_MAX, INT_MAX);
  m.emplace(INT_MAX, 0);

  BOOST_TEST(m.size() == 4);
  BOOST_TEST(m.contains(2));
  BOOST_TEST(m.contains(-2));
  BOOST_TEST(m.contains(54321));
  BOOST_TEST(m.contains(INT_MAX));
  BOOST_TEST(!m.contains(0));
  BOOST_TEST(!m.contains(1));
  BOOST_TEST(!m.contains(-1));
  BOOST_TEST(!m.contains(INT_MIN));
}

BOOST_AUTO_TEST_CASE(flatmap_string_string) {
  dfe::FlatMap<std::string, std::string> m;

  m.emplace("abc", "xyz");
  m.emplace("xyz", "abc");
  m.emplace("xyz", "012");
  m.emplace("012", "asdfsdfgjkl;asdlkjfaiew1-1209vzx-0asdf-as0d8541");

  BOOST_TEST(m.size() == 3);
  BOOST_TEST(m.contains("abc"));
  BOOST_TEST(m.contains("xyz"));
  BOOST_TEST(m.contains("012"));
  BOOST_TEST(!m.contains(""));
  BOOST_TEST(!m.contains("asdfasdgasdgadgs"));
  BOOST_TEST(!m.contains("0124"));
}
