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
#include <iomanip> 
#include <numbers>
#include <utility>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(ToRadian_WrappingCoverage) {
  using std::numbers::pi;

  // {input in degrees, expected radians (wrapped to [-pi, pi))}
  std::vector<std::pair<double, double>> testCases = {
      // In-range
      {0.0, 0.0},
      {90.0, pi / 2.0},
      {180.0, -pi},
      {270.0, -pi / 2.0},
      {359.0, -1.0 * UnitConstants::degree},
      {1.0, 1.0 * UnitConstants::degree},
      {45.0, pi / 4.0},
      {135.0, 3.0 * pi / 4.0},
      {225.0, -3.0 * pi / 4.0},
      {315.0, -pi / 4.0},

      // Just below 360
      {360.0 - 1e-10, -1e-10 * UnitConstants::degree},

      // Negative values
      {-1.0, -1.0 * UnitConstants::degree},
      {-90.0, -pi / 2.0},
      {-180.0, -pi},
      {-270.0, pi / 2.0},
      {-360.0, 0.0},
      {-450.0, -pi / 2.0},

      // > 360 values
      {361.0, 1.0 * UnitConstants::degree},
      {450.0, pi / 2.0},
      {720.0, 0.0},
      {810.0, pi / 2.0},
      {1080.0 + 45.0, pi / 4.0},

      // Large positive and negative inputs
      {1260.0, -pi},         
      {-1260.0, -pi},        
      
  };

  for (const auto& [degInput, expectedWrappedRad] : testCases) {
    const double actual = TGeoSurfaceConverter::toRadian(degInput);
    CHECK_CLOSE_ABS(actual, expectedWrappedRad, s_epsilon);
  }
}

}  // namespace Acts::Test
