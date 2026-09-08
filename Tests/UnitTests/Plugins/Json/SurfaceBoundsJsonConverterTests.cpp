// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsPlugins/Json/SurfaceBoundsJsonConverter.hpp"

#include <algorithm>
#include <fstream>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(SurfaceBoundsRoundTripTests) {
  std::ofstream out;

  // As all SurfaceBounds have the same streaming API only a one is
  // tested here, all others are tests are identical

  auto rectangeRef = std::make_shared<const RectangleBounds>(4., 6.);
  nlohmann::json rectangleOut =
      SurfaceBoundsJsonConverter::toJson(*rectangeRef);
  out.open("RectangleBounds.json");
  out << rectangleOut.dump(2);
  out.close();

  auto in = std::ifstream("RectangleBounds.json",
                          std::ifstream::in | std::ifstream::binary);
  BOOST_CHECK(in.good());
  nlohmann::json rectangleIn;
  in >> rectangleIn;
  in.close();

  auto rectangleTest =
      SurfaceBoundsJsonConverter::fromJson<RectangleBounds>(rectangleIn);

  BOOST_CHECK(rectangeRef->values() == rectangleTest->values());
}

BOOST_AUTO_TEST_CASE(SurfaceBoundsValueCountTests) {
  // CylinderBounds gained eBevelMinZ/eBevelMaxZ after this format was already
  // in use, so payloads written before that carry four values, not six.
  nlohmann::json jShort;
  jShort["type"] = "CylinderBounds";
  jShort["values"] = std::vector<double>{30., 200., std::numbers::pi, 0.};
  BOOST_CHECK_THROW(
      SurfaceBoundsJsonConverter::fromJson<CylinderBounds>(jShort),
      std::invalid_argument);

  nlohmann::json jLong;
  jLong["type"] = "RectangleBounds";
  jLong["values"] = std::vector<double>{-4., -6., 4., 6., 0.};
  BOOST_CHECK_THROW(
      SurfaceBoundsJsonConverter::fromJson<RectangleBounds>(jLong),
      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
