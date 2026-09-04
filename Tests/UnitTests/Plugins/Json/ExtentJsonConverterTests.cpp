// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsPlugins/Json/ExtentJsonConverter.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(ExtentRoundtripTests) {
  Extent e;
  e.set(AxisDirection::AxisR, 0, 200);
  e.set(AxisDirection::AxisZ, -50, 50);

  nlohmann::json j;
  j["extent"] = e;

  std::cout << j.dump(2) << std::endl;

  Extent eIn = j["extent"];

  CHECK_CLOSE_ABS(eIn.min(AxisDirection::AxisR), e.min(AxisDirection::AxisR),
                  10e-5);
  CHECK_CLOSE_ABS(eIn.max(AxisDirection::AxisR), e.max(AxisDirection::AxisR),
                  10e-5);
  CHECK_CLOSE_ABS(eIn.min(AxisDirection::AxisZ), e.min(AxisDirection::AxisZ),
                  10e-5);
  CHECK_CLOSE_ABS(eIn.max(AxisDirection::AxisZ), e.max(AxisDirection::AxisZ),
                  10e-5);
}

BOOST_AUTO_TEST_CASE(ExtentEnvelopeRoundtripTests) {
  ExtentEnvelope env;
  env[AxisDirection::AxisR] = {1., 2.};
  env[AxisDirection::AxisZ] = {3., 4.};

  nlohmann::json j;
  j["envelope"] = env;

  ExtentEnvelope envIn = j["envelope"];

  BOOST_CHECK_EQUAL(envIn[AxisDirection::AxisR][0],
                    env[AxisDirection::AxisR][0]);
  BOOST_CHECK_EQUAL(envIn[AxisDirection::AxisR][1],
                    env[AxisDirection::AxisR][1]);
  BOOST_CHECK_EQUAL(envIn[AxisDirection::AxisZ][0],
                    env[AxisDirection::AxisZ][0]);
  BOOST_CHECK_EQUAL(envIn[AxisDirection::AxisZ][1],
                    env[AxisDirection::AxisZ][1]);
  // Axes not set should remain zero
  BOOST_CHECK_EQUAL(envIn[AxisDirection::AxisX][0], 0.);
  BOOST_CHECK_EQUAL(envIn[AxisDirection::AxisX][1], 0.);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
