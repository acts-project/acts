// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "ActsFatras/Digitization/Segmentation.hpp"

BOOST_AUTO_TEST_SUITE(Digitization)

BOOST_AUTO_TEST_CASE(CartesianSegmentationTests) {
  using namespace ActsFatras;
  using namespace Acts;

  CartesianSegmentation segmentation(
      ProtoAxis(AxisBoundaryType::Bound, -10, 10, 20),
      ProtoAxis(AxisBoundaryType::Bound, -12, 12, 24));

  // Test binning - 0 is underflow bin
  auto bin = segmentation.bin(Vector2(0.1, 0.1));
  BOOST_CHECK_EQUAL(bin[0], 10);
  BOOST_CHECK_EQUAL(bin[1], 12);
  // Test position
  auto pos = segmentation.position(bin);
  CHECK_CLOSE_ABS(pos.x(), 0.5, 1e-5);
  CHECK_CLOSE_ABS(pos.y(), 0.5, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()
