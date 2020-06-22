// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <iterator>
#include <stdexcept>

#include "Acts/Geometry/HierarchicalGeometryContainer.hpp"

namespace {

using Acts::GeometryID;

struct Thing {
  GeometryID id;
  double value = 1.0;

  Thing(GeometryID i, double v) : id(i), value(v) {}

  constexpr auto geometryId() const { return id; }
};

using Container = Acts::HierarchicalGeometryContainer<Thing>;

// helper function to create geometry ids
GeometryID makeId(int volume, int layer = 0, int sensitive = 0) {
  return GeometryID().setVolume(volume).setLayer(layer).setSensitive(sensitive);
}

}  // namespace

// check that an entry exists for the query and compare the id
#define CHECK_ENTRY(container, query, compare) \
  do {                                         \
    auto ret = container.find(query);          \
    BOOST_TEST(ret != container.end());        \
    if (ret != container.end()) {              \
      BOOST_TEST(ret->id == compare);          \
    }                                          \
  } while (false)

BOOST_TEST_DONT_PRINT_LOG_VALUE(Container::Iterator)
BOOST_TEST_DONT_PRINT_LOG_VALUE(Thing)

BOOST_AUTO_TEST_SUITE(HierarchicalGeometryContainer)

BOOST_AUTO_TEST_CASE(ConstructDefault) {
  Container c;
  BOOST_TEST(c.begin() == c.end());
  BOOST_TEST(c.empty());
  BOOST_TEST(c.size() == 0u);
}

BOOST_AUTO_TEST_CASE(ConstructNonUnique) {
  std::vector<Thing> elements = {
      {makeId(2, 4, 6), 1.0},
      {makeId(3, 5), 1.0},
      {makeId(3), 1.0},
      // duplicate identifier
      {makeId(2, 4, 6), 2.0},
  };
  BOOST_CHECK_THROW(Container(std::move(elements)), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ConstructInitializerList) {
  Container c = {
      {makeId(0, 1, 2), 1.0},
      {makeId(3, 4), 3.0},
      {makeId(2), 2.0},
  };
  BOOST_TEST(std::next(c.begin(), 3) == c.end());
  BOOST_TEST(not c.empty());
  BOOST_TEST(c.size() == 3u);
  // only test that all elements are there; failure test are below
  CHECK_ENTRY(c, makeId(0, 1, 2), makeId(0, 1, 2));
  CHECK_ENTRY(c, makeId(2), makeId(2));
  CHECK_ENTRY(c, makeId(3, 4), makeId(3, 4));
}

BOOST_AUTO_TEST_CASE(IndexBasedAccess) {
  Container c({
      {makeId(1, 2, 3), 1.0},
      {makeId(3, 4, 5), 2.0},
      {makeId(3, 5), 3.0},
      {makeId(4, 5, 7), 4.0},
  });

  BOOST_TEST(not c.empty());
  BOOST_TEST(c.size() == 4u);
  // this test both that the index-based access works and that the identifier
  // stored in the container matches the one from the stored element
  // NOTE order is undefined and should not be tested
  for (auto i = c.size(); 0 < i--;) {
    BOOST_TEST(c.idAt(i) == c.valueAt(i).id);
  }
  // test that invalid inputs actually fail
  BOOST_CHECK_THROW(c.idAt(c.size()), std::out_of_range);
  BOOST_CHECK_THROW(c.valueAt(c.size()), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(Find) {
  Container c({
      // entry for volume 2, layer 4, sensitive 6
      {makeId(2, 4, 6), -23.0},
      // entry for volume 2, layer 8
      {makeId(2, 8), 5.0},
      // entry for the whole volume 2
      {makeId(2), 1.0},
      // entry for volume 12, layer 16
      // NOTE no entry for the volume as a whole
      {makeId(12, 16), -1.0},
  });

  // basic checks
  BOOST_TEST(std::next(c.begin(), 4u) == c.end());
  BOOST_TEST(not c.empty());
  BOOST_TEST(c.size() == 4u);

  // find existing sensitive
  CHECK_ENTRY(c, makeId(2, 4, 6), makeId(2, 4, 6));
  // find existing layer
  CHECK_ENTRY(c, makeId(2, 8), makeId(2, 8));
  // find existing volume
  CHECK_ENTRY(c, makeId(2), makeId(2));
  // find existing layer
  CHECK_ENTRY(c, makeId(12, 16), makeId(12, 16));

  // find non-existing sensitive, which has a set volume
  CHECK_ENTRY(c, makeId(2, 4, 7), makeId(2));
  // find non-existing layer id, which has a set volume
  CHECK_ENTRY(c, makeId(2, 13), makeId(2));
  // find non-existing sensitive id, which has a set layer
  CHECK_ENTRY(c, makeId(2, 8, 13), makeId(2, 8));
  // find non-existing sensitive, which has a set layer
  CHECK_ENTRY(c, makeId(12, 16, 20), makeId(12, 16));

  // find non-existing sensitive, which has no higher hierarchy set
  BOOST_TEST(c.find(makeId(3, 5, 7)) == c.end());
  // find non-existing layer, which has no higher hierarchy set
  BOOST_TEST(c.find(makeId(3, 5)) == c.end());
  // find non-existing volume
  BOOST_TEST(c.find(makeId(3)) == c.end());
}

BOOST_AUTO_TEST_SUITE_END()
