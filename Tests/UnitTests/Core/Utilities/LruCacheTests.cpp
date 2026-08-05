// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/detail/LruCache.hpp"

#include <string>

namespace {
using Cache = Acts::detail::LruCache<std::string, int>;
}

BOOST_AUTO_TEST_SUITE(LruCacheTests)

BOOST_AUTO_TEST_CASE(ZeroCapacityThrows) {
  BOOST_CHECK_THROW(Cache{0}, std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MissReturnsNullopt) {
  Cache c{2};
  BOOST_CHECK(!c.get("x").has_value());
}

BOOST_AUTO_TEST_CASE(PutAndGet) {
  Cache c{2};
  c.put("a", 1);
  auto v = c.get("a");
  BOOST_REQUIRE(v.has_value());
  BOOST_CHECK_EQUAL(*v, 1);
  BOOST_CHECK_EQUAL(c.size(), 1u);
}

BOOST_AUTO_TEST_CASE(UpdateExistingKey) {
  Cache c{2};
  c.put("a", 1);
  c.put("a", 99);
  BOOST_CHECK_EQUAL(*c.get("a"), 99);
  BOOST_CHECK_EQUAL(c.size(), 1u);
}

BOOST_AUTO_TEST_CASE(EvictsLRU) {
  Cache c{2};
  c.put("a", 1);
  c.put("b", 2);
  c.put("c", 3);  // evicts "a" (LRU)
  BOOST_CHECK(!c.get("a").has_value());
  BOOST_CHECK_EQUAL(*c.get("b"), 2);
  BOOST_CHECK_EQUAL(*c.get("c"), 3);
  BOOST_CHECK_EQUAL(c.size(), 2u);
}

BOOST_AUTO_TEST_CASE(GetPromotesToMRU) {
  Cache c{2};
  c.put("a", 1);
  c.put("b", 2);
  c.get("a");     // promotes "a" to MRU
  c.put("c", 3);  // evicts "b" (now LRU), not "a"
  BOOST_CHECK_EQUAL(*c.get("a"), 1);
  BOOST_CHECK(!c.get("b").has_value());
  BOOST_CHECK_EQUAL(*c.get("c"), 3);
}

BOOST_AUTO_TEST_CASE(CapacityOne) {
  Cache c{1};
  c.put("a", 1);
  c.put("b", 2);
  BOOST_CHECK(!c.get("a").has_value());
  BOOST_CHECK_EQUAL(*c.get("b"), 2);
  BOOST_CHECK_EQUAL(c.size(), 1u);
}

BOOST_AUTO_TEST_CASE(SizeAndCapacity) {
  Cache c{3};
  BOOST_CHECK_EQUAL(c.capacity(), 3u);
  BOOST_CHECK_EQUAL(c.size(), 0u);
  c.put("a", 1);
  BOOST_CHECK_EQUAL(c.size(), 1u);
  c.put("b", 2);
  c.put("c", 3);
  BOOST_CHECK_EQUAL(c.size(), 3u);
  c.put("d", 4);  // evicts one
  BOOST_CHECK_EQUAL(c.size(), 3u);
}

BOOST_AUTO_TEST_SUITE_END()
