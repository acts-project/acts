// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/detail/UtilityFunctions.hpp"

#include <iterator>
#include <list>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(insert_vector) {
  std::vector<std::size_t> coll;
  Acts::detail::insert(coll, 2ul);
  BOOST_CHECK(coll.size() == 1ul);
  Acts::detail::insert(coll, 5ul);
  BOOST_CHECK(coll.size() == 2ul);
  std::size_t val = 1ul;
  Acts::detail::insert(coll, val);
  BOOST_CHECK(coll.size() == 3ul);

  BOOST_CHECK_EQUAL(coll[0], 2ul);
  BOOST_CHECK_EQUAL(coll[1], 5ul);
  BOOST_CHECK_EQUAL(coll[2], 1ul);
}

BOOST_AUTO_TEST_CASE(insert_list) {
  std::list<std::size_t> coll;
  Acts::detail::insert(coll, 2ul);
  BOOST_CHECK(coll.size() == 1ul);
  Acts::detail::insert(coll, 5ul);
  BOOST_CHECK(coll.size() == 2ul);
  std::size_t val = 1ul;
  Acts::detail::insert(coll, val);
  BOOST_CHECK(coll.size() == 3ul);

  BOOST_CHECK_EQUAL(coll.front(), 2ul);
  coll.pop_front();
  BOOST_CHECK_EQUAL(coll.front(), 5ul);
  coll.pop_front();
  BOOST_CHECK_EQUAL(coll.front(), 1ul);
  coll.pop_front();
}

BOOST_AUTO_TEST_CASE(insert_set) {
  std::set<std::size_t> coll;
  Acts::detail::insert(coll, 2ul);
  BOOST_CHECK(coll.size() == 1ul);
  Acts::detail::insert(coll, 5ul);
  BOOST_CHECK(coll.size() == 2ul);
  std::size_t val = 1ul;
  Acts::detail::insert(coll, val);
  BOOST_CHECK(coll.size() == 3ul);

  BOOST_CHECK(coll.find(2ul) != coll.end());
  BOOST_CHECK(coll.find(5ul) != coll.end());
  BOOST_CHECK(coll.find(1ul) != coll.end());
}

BOOST_AUTO_TEST_CASE(insert_unordered_set) {
  std::unordered_set<std::size_t> coll;
  Acts::detail::insert(coll, 2ul);
  BOOST_CHECK(coll.size() == 1ul);
  Acts::detail::insert(coll, 5ul);
  BOOST_CHECK(coll.size() == 2ul);
  std::size_t val = 1ul;
  Acts::detail::insert(coll, val);
  BOOST_CHECK(coll.size() == 3ul);

  BOOST_CHECK(coll.find(2ul) != coll.end());
  BOOST_CHECK(coll.find(5ul) != coll.end());
  BOOST_CHECK(coll.find(1ul) != coll.end());
}

BOOST_AUTO_TEST_CASE(insert_back_iterator) {
  std::vector<std::size_t> coll;
  std::back_insert_iterator<decltype(coll)> backItr(coll);
  Acts::detail::insert(backItr, 2ul);
  BOOST_CHECK(coll.size() == 1ul);
  Acts::detail::insert(backItr, 5ul);
  BOOST_CHECK(coll.size() == 2ul);
  std::size_t val = 1ul;
  Acts::detail::insert(backItr, val);
  BOOST_CHECK(coll.size() == 3ul);

  BOOST_CHECK_EQUAL(coll[0], 2ul);
  BOOST_CHECK_EQUAL(coll[1], 5ul);
  BOOST_CHECK_EQUAL(coll[2], 1ul);
}

BOOST_AUTO_TEST_CASE(emplace_vector) {
  struct A {
    A(std::size_t a, float b, std::string c)
        : valA(a), valB(b), valC(std::move(c)) {}

    std::size_t valA;
    float valB;
    std::string valC;
  };

  std::vector<A> coll;
  Acts::detail::emplace(coll, 2ul, 4.5f, "A");
  BOOST_CHECK(coll.size() == 1ul);
  std::size_t a = 5ul;
  float b = 0.25f;
  std::string c = "B";
  Acts::detail::emplace(coll, a, b, c);
  BOOST_CHECK(coll.size() == 2ul);

  BOOST_CHECK_EQUAL(coll[0].valA, 2ul);
  BOOST_CHECK_EQUAL(coll[0].valB, 4.5f);
  BOOST_CHECK_EQUAL(coll[0].valC, "A");

  BOOST_CHECK_EQUAL(coll[1].valA, a);
  BOOST_CHECK_EQUAL(coll[1].valB, b);
  BOOST_CHECK_EQUAL(coll[1].valC, c);
}

}  // namespace Acts::Test
