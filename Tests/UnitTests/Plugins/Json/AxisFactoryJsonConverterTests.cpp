// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisFactory.hpp"
#include "Acts/Utilities/MultiAxisFactory.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/AxisFactoryJsonConverter.hpp"

#include <stdexcept>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(AxisFactoryJsonRoundTrip) {
  using enum AxisBoundaryType;
  using enum AxisDirection;

  std::vector<AxisFactory> descriptions = {
      AxisFactory::Equidistant(Bound, 0., 10., 5),
      AxisFactory::Equidistant(Closed, -3., 3., 12, AxisPhi),
      AxisFactory::Variable(Open, {0., 1., 4., 10.}),
      AxisFactory::Variable(Bound, {-1., 0., 2.}, AxisZ),
      AxisFactory::DeferredEquidistant(20),
      AxisFactory::DeferredEquidistant(8, AxisRPhi),
      AxisFactory::DeferredVariable({0., 0.1, 0.5, 1.}),
      AxisFactory::DeferredVariable({0., 0.25, 1.}, AxisR)};

  for (const AxisFactory& axisFactory : descriptions) {
    nlohmann::json j = AxisFactoryJsonConverter::toJson(axisFactory);
    AxisFactory read = AxisFactoryJsonConverter::fromJson(j);
    BOOST_CHECK(read == axisFactory);
    // The direction key is only written when a direction is set
    BOOST_CHECK_EQUAL(j.contains("direction"),
                      axisFactory.direction().has_value());
  }

  // Unknown type tag is rejected
  nlohmann::json jInvalid = {{"type", "unknown"}};
  BOOST_CHECK_THROW(AxisFactoryJsonConverter::fromJson(jInvalid),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MultiAxisFactoryJsonRoundTrip) {
  using enum AxisBoundaryType;
  using enum AxisDirection;

  MultiAxisFactory oneD({AxisFactory::DeferredEquidistant(8, AxisZ)});
  nlohmann::json j1 = MultiAxisFactoryJsonConverter::toJson(oneD);
  BOOST_CHECK(MultiAxisFactoryJsonConverter::fromJson(j1) == oneD);

  MultiAxisFactory twoD({AxisFactory::DeferredEquidistant(4, AxisRPhi),
                         AxisFactory::DeferredVariable({0., 0.5, 1.}, AxisZ)});
  nlohmann::json j2 = MultiAxisFactoryJsonConverter::toJson(twoD);
  BOOST_CHECK_EQUAL(j2.size(), 2u);
  BOOST_CHECK(MultiAxisFactoryJsonConverter::fromJson(j2) == twoD);

  // An empty axis list is rejected
  BOOST_CHECK_THROW(
      MultiAxisFactoryJsonConverter::fromJson(nlohmann::json::array()),
      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
