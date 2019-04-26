// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE TrackingGeometry Tests

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"

namespace Acts {
namespace Test {

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  BOOST_AUTO_TEST_CASE(CylindricalTrackingGeometryTest)
  {
    CylindricalTrackingGeometry cGeometry(tgContext);
    auto                        tGeometry = cGeometry();
    BOOST_CHECK_NE(tGeometry, nullptr);
  }

  BOOST_AUTO_TEST_CASE(CubicTrackingGeometryTest)
  {
    CubicTrackingGeometry cGeometry(tgContext);
    auto                  tGeometry = cGeometry();
    BOOST_CHECK_NE(tGeometry, nullptr);
  }

}  // namespace Acts
}  // namespace Test
