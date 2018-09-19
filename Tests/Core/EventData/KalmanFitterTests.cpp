// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include "Acts/Detector/TrackingGeometry.hpp"
#include "DetectorBuild.hpp"

#include <vector>

namespace Acts {
namespace Test {

  ///
  /// @brief Unit test for Kalman fitter
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_initialization)
  {
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

    BOOST_TEST(true);
  }

}  // namespace Test
}  // namespace Acts
