// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/AxisSpecJsonConverter.hpp"

#include <stdexcept>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(AxisSpecJsonRoundTrip) {
  using enum AxisBoundaryType;
  using enum AxisDirection;

  std::vector<AxisSpec> specs = {
      AxisSpec::Equidistant(5, 0., 10., Bound),
      AxisSpec::Equidistant(12, -3., 3., Closed, AxisPhi),
      AxisSpec::Variable({0., 1., 4., 10.}, Open),
      AxisSpec::Variable({-1., 0., 2.}, Bound, AxisZ),
      AxisSpec::DeferredEquidistant(20),
      AxisSpec::DeferredEquidistant(8, AxisRPhi),
      AxisSpec::DeferredVariable({0., 0.1, 0.5, 1.}),
      AxisSpec::DeferredVariable({0., 0.25, 1.}, std::nullopt, AxisR),
      // Partially specified: a range without a boundary type and the other
      // way round
      AxisSpec::Equidistant(6, 0., 1.),
      AxisSpec::Equidistant(6, std::nullopt, std::nullopt, Bound),
      AxisSpec::Variable({0., 1., 2.}),
      AxisSpec::DeferredVariable({0., 0.5, 1.}, Closed)};

  for (const AxisSpec& axisSpec : specs) {
    nlohmann::json j = AxisSpecJsonConverter::toJson(axisSpec);
    AxisSpec read = AxisSpecJsonConverter::fromJson(j);
    BOOST_CHECK(read == axisSpec);
    // The direction key is only written when a direction is set
    BOOST_CHECK_EQUAL(j.contains("direction"),
                      axisSpec.direction().has_value());
  }

  // Unknown type tag is rejected
  nlohmann::json jInvalid = {{"type", "unknown"}};
  BOOST_CHECK_THROW(AxisSpecJsonConverter::fromJson(jInvalid),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MultiAxisSpecJsonRoundTrip) {
  using enum AxisBoundaryType;
  using enum AxisDirection;

  MultiAxisSpec oneD({AxisSpec::DeferredEquidistant(8, AxisZ)});
  nlohmann::json j1 = MultiAxisSpecJsonConverter::toJson(oneD);
  BOOST_CHECK(MultiAxisSpecJsonConverter::fromJson(j1) == oneD);

  MultiAxisSpec twoD(
      {AxisSpec::DeferredEquidistant(4, AxisRPhi),
       AxisSpec::DeferredVariable({0., 0.5, 1.}, std::nullopt, AxisZ)});
  nlohmann::json j2 = MultiAxisSpecJsonConverter::toJson(twoD);
  BOOST_CHECK_EQUAL(j2.size(), 2u);
  BOOST_CHECK(MultiAxisSpecJsonConverter::fromJson(j2) == twoD);

  // An empty axis list is rejected
  BOOST_CHECK_THROW(
      MultiAxisSpecJsonConverter::fromJson(nlohmann::json::array()),
      std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
