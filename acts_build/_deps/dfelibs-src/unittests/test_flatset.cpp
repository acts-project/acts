// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for dfe::FlatSet

#include <boost/test/unit_test.hpp>

#include "dfe/dfe_flat.hpp"

using dfe::FlatSet;

BOOST_AUTO_TEST_CASE(flatset_int) {
  FlatSet<int> set;

  set.insert_or_assign(1);
  set.insert_or_assign(-3);
  set.insert_or_assign(12);
  set.insert_or_assign(4);

  BOOST_CHECK(set.size() == 4);
  BOOST_CHECK(set.contains(-3));
  BOOST_CHECK(set.contains(1));
  BOOST_CHECK(set.contains(4));
  BOOST_CHECK(set.contains(12));
  BOOST_CHECK(!set.contains(-5));
  BOOST_CHECK(!set.contains(23));
  BOOST_CHECK(!set.contains(42));

  // test checked access
  BOOST_CHECK_NO_THROW(set.at(4));
  BOOST_CHECK_NO_THROW(set.at(12));
  BOOST_CHECK_THROW(set.at(-5), std::out_of_range);
  BOOST_CHECK_THROW(set.at(231), std::out_of_range);

  // test iteration and ordering
  auto it = set.begin();
  BOOST_CHECK(*it == -3);
  ++it;
  BOOST_CHECK(*it == 1);
  ++it;
  BOOST_CHECK(*it == 4);
  ++it;
  BOOST_CHECK(*it == 12);
  ++it;
  BOOST_CHECK(it == set.end());

  // try to re-insert existing objects
  for (const auto& x : set) {
    auto size_before = set.size();
    set.insert_or_assign(x);
    BOOST_CHECK(set.size() == size_before);
  }
}

struct Thing {
  int index;
  float value;
};
struct ThingComparator {
  bool operator()(const Thing& a, const Thing& b) { return a.index < b.index; }
  bool operator()(int aindex, const Thing& b) { return aindex < b.index; }
  bool operator()(const Thing& a, int bindex) { return a.index < bindex; }
};

BOOST_AUTO_TEST_CASE(flatset_custom_compare) {
  FlatSet<Thing, ThingComparator> set;

  set.insert_or_assign({12, 0.25});
  set.insert_or_assign({23, 1.25});
  BOOST_CHECK(set.size() == 2);

  // check replacement of objects w/ same identity
  BOOST_CHECK(set.find(23) != set.end());
  BOOST_CHECK(set.find(23)->index == 23);
  BOOST_CHECK(set.find(23)->value == 1.25);
  set.insert_or_assign({23, 5000.5});
  BOOST_CHECK(set.size() == 2);
  BOOST_CHECK(set.find(23) != set.end());
  BOOST_CHECK(set.find(23)->index == 23);
  BOOST_CHECK(set.find(23)->value == 5000.5);

  // identity is defined by the compare functionality, not by value
  BOOST_CHECK(set.contains(12));
  BOOST_CHECK(set.contains(Thing{12, 0.0}));
  BOOST_CHECK(set.contains(Thing{12, -1.23}));
  BOOST_CHECK(set.contains(23));
  BOOST_CHECK(set.contains(Thing{23, 0.0}));
  BOOST_CHECK(set.contains(Thing{23, -1.23}));
  BOOST_CHECK(!set.contains(4));
  BOOST_CHECK(!set.contains(Thing{4, 1.45}));
  BOOST_CHECK(!set.contains(Thing{27, -1.23}));
}
