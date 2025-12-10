// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"

using namespace Acts;
using namespace ActsExamples;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(AlignmentSuite)

BOOST_AUTO_TEST_CASE(Single_sensitive_Test) {
  auto geoId = GeometryIdentifier().withSensitive(1);
  AlignmentGroup aGroup("test", {geoId});
  BOOST_CHECK(aGroup.has(geoId));

  auto badId = GeometryIdentifier().withSensitive(2);
  BOOST_CHECK(!aGroup.has(badId));
}

BOOST_AUTO_TEST_CASE(Single_layer_Test) {
  auto geoId = GeometryIdentifier().withLayer(1);
  AlignmentGroup aGroup("test", {geoId});
  BOOST_CHECK(aGroup.has(geoId));

  auto badId = GeometryIdentifier().withLayer(2);
  BOOST_CHECK(!aGroup.has(badId));
}

BOOST_AUTO_TEST_CASE(Single_volume_Test) {
  auto geoId = GeometryIdentifier().withVolume(1);
  AlignmentGroup aGroup("test", {geoId});
  BOOST_CHECK(aGroup.has(geoId));

  auto badId = GeometryIdentifier().withVolume(2);
  BOOST_CHECK(!aGroup.has(badId));
}

BOOST_AUTO_TEST_CASE(Hierarchy_test) {
  auto geoId = GeometryIdentifier().withVolume(1);
  AlignmentGroup aGroup("test", {geoId});

  auto geoId2 = GeometryIdentifier().withVolume(1).withLayer(2);
  BOOST_CHECK(aGroup.has(geoId2));

  auto geoId3 = GeometryIdentifier().withVolume(1).withLayer(3);
  BOOST_CHECK(aGroup.has(geoId3));
}

BOOST_AUTO_TEST_CASE(Hierarchy_test_2) {
  auto geoId = GeometryIdentifier().withVolume(1).withLayer(1);
  AlignmentGroup aGroup("test", {geoId});

  auto badId = GeometryIdentifier().withVolume(1);
  BOOST_CHECK(!aGroup.has(badId));
}

BOOST_AUTO_TEST_CASE(Multiple_test) {
  auto geoId1 = GeometryIdentifier().withVolume(1).withLayer(1);
  auto geoId2 = GeometryIdentifier().withVolume(2).withLayer(2);
  auto geoId3 = GeometryIdentifier().withVolume(3).withLayer(3);
  auto geoId4 = GeometryIdentifier().withVolume(4).withLayer(4);

  AlignmentGroup aGroup("test", {geoId1, geoId2, geoId3});

  BOOST_CHECK(aGroup.has(geoId1));
  BOOST_CHECK(aGroup.has(geoId2));
  BOOST_CHECK(aGroup.has(geoId3));
  BOOST_CHECK(!aGroup.has(geoId4));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
