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

BOOST_TEST_DONT_PRINT_LOG_VALUE(Container::Iterator)

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
  {
    auto ret = c.find(makeId(2, 4, 6));
    BOOST_TEST(ret != c.end());
    BOOST_TEST(ret->id == makeId(2, 4, 6));
  }
  // find non-existing sensitive, which has a set volume
  {
    auto ret = c.find(makeId(2, 4, 7));
    BOOST_TEST(ret != c.end());
    BOOST_TEST(ret->id == makeId(2));
  }
  // find non-existing layer id, which has a set volume
  {
    auto ret = c.find(makeId(2, 13));
    BOOST_TEST(ret != c.end());
    BOOST_TEST(ret->id == makeId(2));
  }
  // find non-existing sensitive id, which has a set layer
  {
    auto ret = c.find(makeId(2, 8, 13));
    BOOST_TEST(ret != c.end());
    BOOST_TEST(ret->id == makeId(2, 8));
  }
  // find existing layer
  {
    auto ret = c.find(makeId(2, 8));
    BOOST_TEST(ret != c.end());
    BOOST_TEST(ret->id == makeId(2, 8));
  }
  // find existing volume
  {
    auto ret = c.find(makeId(2));
    BOOST_TEST(ret != c.end());
    BOOST_TEST(ret->id == makeId(2));
  }

  // find non-existing sensitive, which has a set layer
  {
    auto ret = c.find(makeId(12, 16, 20));
    BOOST_TEST(ret != c.end());
    BOOST_TEST(ret->id == makeId(12, 16));
  }
  // find existing layer
  {
    auto ret = c.find(makeId(12, 16));
    BOOST_TEST(ret != c.end());
    BOOST_TEST(ret->id == makeId(12, 16));
  }
  // find non-existing volume, which has only lower hierarchy elements
  BOOST_TEST(c.find(makeId(12)) == c.end());

  // find non-existing sensitive, which has no higher hierarchy set
  BOOST_TEST(c.find(makeId(3, 5, 7)) == c.end());
  // find non-existing layer, which has no higher hierarchy set
  BOOST_TEST(c.find(makeId(3, 5)) == c.end());
  // find non-existing volume
  BOOST_TEST(c.find(makeId(3)) == c.end());
}

BOOST_AUTO_TEST_SUITE_END()
