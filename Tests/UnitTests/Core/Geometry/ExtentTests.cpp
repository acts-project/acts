// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

// Helper
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

// The class to test
#include "Acts/Geometry/Extent.hpp"

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

using namespace UnitLiterals;

namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

/// Unit tests for Polyderon construction & operator +=
BOOST_AUTO_TEST_CASE(ExtentTest) {
  std::vector<Vector3D> vertices = {
      Vector3D(15_mm, -3_mm, -10_mm), Vector3D(18_mm, 0_mm, -10_mm),
      Vector3D(15_mm, 3_mm, -10_mm),  Vector3D(15_mm, -3_mm, 10_mm),
      Vector3D(18_mm, 0_mm, 10_mm),   Vector3D(15_mm, 3_mm, 10_mm)};

  // Create an Extent
  Extent gExt;
  for (const auto& v : vertices) {
    gExt.check(v);
  }

  double phiMin = std::atan2(-3_mm, 15_mm);
  double phiMax = std::atan2(3_mm, 15_mm);
  double rMin = std::sqrt(15_mm * 15_mm + 3_mm * 3_mm);

  CHECK_CLOSE_ABS(gExt.min(binX), 15_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binX), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.min(binY), -3_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binY), 3_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.min(binZ), -10_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binZ), 10_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.min(binR), rMin, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binR), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.min(binPhi), phiMin, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binPhi), phiMax, 1e-6);

  // Create a second Extent
  Extent otherExt;
  otherExt.extend(gExt);

  CHECK_CLOSE_ABS(otherExt.min(binX), 15_mm, 1e-6);
  CHECK_CLOSE_ABS(otherExt.max(binX), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(otherExt.min(binY), -3_mm, 1e-6);
  CHECK_CLOSE_ABS(otherExt.max(binY), 3_mm, 1e-6);
  CHECK_CLOSE_ABS(otherExt.min(binZ), -10_mm, 1e-6);
  CHECK_CLOSE_ABS(otherExt.max(binZ), 10_mm, 1e-6);
  CHECK_CLOSE_ABS(otherExt.min(binR), rMin, 1e-6);
  CHECK_CLOSE_ABS(otherExt.max(binR), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(otherExt.min(binPhi), phiMin, 1e-6);
  CHECK_CLOSE_ABS(otherExt.max(binPhi), phiMax, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts