// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"

#include "Acts/Utilities/BinAdjustment.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <cmath>

namespace Acts {
namespace Test {

namespace tt = boost::test_tools;

// Test Radial
BOOST_AUTO_TEST_CASE(BinAdjustment_Radial) {
  RadialBounds bound(50, 75, M_PI, 0);
  BinUtility bu;
  bu += BinUtility(1, 0, 1, Acts::open, Acts::binR);
  bu += BinUtility(1, 0, 1, Acts::closed, Acts::binPhi);

  BinUtility buAdjust = adjustBinUtility(bu, bound, Transform3D::Identity());

  BOOST_CHECK_EQUAL(buAdjust.binningData()[0].min, 50);
  BOOST_CHECK_EQUAL(buAdjust.binningData()[0].max, 75);
  BOOST_CHECK_EQUAL(buAdjust.binningData()[1].min, float(-M_PI));
  BOOST_CHECK_EQUAL(buAdjust.binningData()[1].max, float(M_PI));
}

// Test Cylinder
BOOST_AUTO_TEST_CASE(BinAdjustment_Cylinder) {
  CylinderBounds bound(25, 50, M_PI / 4, 0);
  BinUtility bu;
  bu += BinUtility(1, 0, 1, Acts::open, Acts::binPhi);
  bu += BinUtility(1, 0, 1, Acts::open, Acts::binZ);

  BinUtility buAdjust = adjustBinUtility(bu, bound, Transform3D::Identity());

  BOOST_CHECK_EQUAL(buAdjust.binningData()[0].min, float(-M_PI / 4));
  BOOST_CHECK_EQUAL(buAdjust.binningData()[0].max, float(M_PI / 4));
  BOOST_CHECK_EQUAL(buAdjust.binningData()[1].min, -50);
  BOOST_CHECK_EQUAL(buAdjust.binningData()[1].max, 50);
}

}  // namespace Test
}  // namespace Acts
