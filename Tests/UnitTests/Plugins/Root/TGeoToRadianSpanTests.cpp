// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/Root/TGeoSurfaceConverter.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(ToRadianSpan_Correctness) {
  using std::numbers::pi;

  struct TestCase {
    double deg1;
    double deg2;
    double expectedHalfSpan;
    double expectedAvg;
  };

  std::vector<TestCase> testCases = {
      // Simple cases
      {0.0, 90.0, pi / 4.0, pi / 4.0},
      {90.0, 180.0, pi / 4.0, 3.0 * pi / 4.0},
      {0.0, 360.0, pi, 0.0},
      {360.0, 720.0, pi, 0.0},
      {-180.0, 180.0, 0.0, pi},

      // Wraparound case (crosses 0°)
      {315.0, 45.0, pi / 4.0, 0.0},

      // Reverse direction (same as above, but in opposite order)
      {45.0, 315.0, 3.0 * pi / 4.0, pi},

      // Small wraparound
      {359.0, 1.0, UnitConstants::degree, 0.0},

      // Large positive range
      {0.0, 1125.0, pi / 8.0, pi / 8.0},

      // Large negative to large positive
      {-720.0, 90.0, pi / 4.0, pi / 4.0},

      // Custom case 
      {316.79, 403.21, 43.21 * UnitConstants::degree * 0.5,
       Acts::detail::radian_sym((316.79 + 0.5 * 43.21) * UnitConstants::degree)},
  };

  for (const auto& [deg1, deg2, expectedHalfSpan, expectedAvg] : testCases) {
    const auto [actualHalfSpan, actualAvg] =
        Acts::TGeoSurfaceConverter::toRadianSpan(deg1, deg2);

    BOOST_TEST_CONTEXT("deg1 = " << deg1 << ", deg2 = " << deg2) {
      CHECK_CLOSE_ABS(actualHalfSpan, expectedHalfSpan, s_epsilon);
      CHECK_CLOSE_ABS(actualAvg, expectedAvg, s_epsilon);
    }
  }
}

}  // namespace Acts::Test
