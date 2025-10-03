// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(ReadingWritingFull) {
  GeometryIdentifier geoId = GeometryIdentifier()
                                 .withVolume(1)
                                 .withBoundary(2)
                                 .withLayer(3)
                                 .withApproach(4)
                                 .withSensitive(5)
                                 .withExtra(6);

  nlohmann::json j = GeometryIdentifierJsonConverter::encodeIdentifier(geoId);

  BOOST_CHECK_EQUAL(j["volume"], 1);
  BOOST_CHECK_EQUAL(j["boundary"], 2);
  BOOST_CHECK_EQUAL(j["layer"], 3);
  BOOST_CHECK_EQUAL(j["approach"], 4);
  BOOST_CHECK_EQUAL(j["sensitive"], 5);
  BOOST_CHECK_EQUAL(j["extra"], 6);

  // or simply
  nlohmann::json j2 = nlohmann::json(geoId);
  BOOST_CHECK_EQUAL(j2["volume"], 1);
  BOOST_CHECK_EQUAL(j2["boundary"], 2);
  BOOST_CHECK_EQUAL(j2["layer"], 3);
  BOOST_CHECK_EQUAL(j2["approach"], 4);
  BOOST_CHECK_EQUAL(j2["sensitive"], 5);
  BOOST_CHECK_EQUAL(j2["extra"], 6);

  GeometryIdentifier geoId2 =
      GeometryIdentifierJsonConverter::decodeIdentifier(j);
  BOOST_CHECK_EQUAL(geoId2.volume(), 1);
  BOOST_CHECK_EQUAL(geoId2.boundary(), 2);
  BOOST_CHECK_EQUAL(geoId2.layer(), 3);
  BOOST_CHECK_EQUAL(geoId2.approach(), 4);
  BOOST_CHECK_EQUAL(geoId2.sensitive(), 5);
  BOOST_CHECK_EQUAL(geoId2.extra(), 6);
  // or simply
  GeometryIdentifier geoId3 = j2.get<GeometryIdentifier>();
  BOOST_CHECK_EQUAL(geoId3.volume(), 1);
  BOOST_CHECK_EQUAL(geoId3.boundary(), 2);
  BOOST_CHECK_EQUAL(geoId3.layer(), 3);
  BOOST_CHECK_EQUAL(geoId3.approach(), 4);
  BOOST_CHECK_EQUAL(geoId3.sensitive(), 5);
  BOOST_CHECK_EQUAL(geoId3.extra(), 6);
}

BOOST_AUTO_TEST_CASE(ReadingWritingCompact) {
  GeometryIdentifier geoId = GeometryIdentifier()
                                 .withVolume(1)
                                 .withBoundary(2)
                                 .withLayer(3)
                                 .withApproach(4)
                                 .withSensitive(5)
                                 .withExtra(6);

  nlohmann::json j =
      GeometryIdentifierJsonConverter::encodeIdentifier(geoId, true);

  auto value = j.get<GeometryIdentifier::Value>();
  BOOST_CHECK_EQUAL(value, geoId.value());

  // Rely on the auto-detection of the compact reading
  GeometryIdentifier geoIdCompact = j.get<GeometryIdentifier>();
  BOOST_CHECK_EQUAL(geoIdCompact.volume(), 1);
  BOOST_CHECK_EQUAL(geoIdCompact.boundary(), 2);
  BOOST_CHECK_EQUAL(geoIdCompact.layer(), 3);
  BOOST_CHECK_EQUAL(geoIdCompact.approach(), 4);
  BOOST_CHECK_EQUAL(geoIdCompact.sensitive(), 5);
  BOOST_CHECK_EQUAL(geoIdCompact.extra(), 6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
