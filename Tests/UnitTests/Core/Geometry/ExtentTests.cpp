// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace Acts {

using namespace UnitLiterals;

namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

/// Unit tests for Polyderon construction & operator +=
BOOST_AUTO_TEST_CASE(ExtentTest) {
  std::vector<Vector3> vertices = {
      Vector3(15_mm, -3_mm, -10_mm), Vector3(18_mm, 0_mm, -10_mm),
      Vector3(15_mm, 3_mm, -10_mm),  Vector3(15_mm, -3_mm, 10_mm),
      Vector3(18_mm, 0_mm, 10_mm),   Vector3(15_mm, 3_mm, 10_mm)};

  // Create an Extent / without envelope
  Extent gExt;
  for (const auto& v : vertices) {
    gExt.extend(v);
  }

  double phiMin = std::atan2(-3_mm, 15_mm);
  double phiMax = std::atan2(3_mm, 15_mm);
  double rMin = std::hypot(15_mm, 3_mm);

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

  // Call with histogram filling
  Extent gExtHist;
  for (const auto& v : vertices) {
    gExtHist.extend(v, {binX}, false, true);
  }
  const auto& vHist = gExtHist.valueHistograms();
  auto xVals = vHist[binX];

  BOOST_CHECK_EQUAL(xVals.size(), 6u);
  std::vector<ActsScalar> reference = {15_mm, 18_mm, 15_mm,
                                       15_mm, 18_mm, 15_mm};
  BOOST_CHECK(xVals == reference);

  // Call with ieterator range
  Extent gExtItr;
  gExtItr.extend(vertices.begin(), vertices.end());
  CHECK_CLOSE_ABS(gExtItr.min(binX), 15_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.max(binX), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.min(binY), -3_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.max(binY), 3_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.min(binZ), -10_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.max(binZ), 10_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.min(binR), rMin, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.max(binR), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.min(binPhi), phiMin, 1e-6);
  CHECK_CLOSE_ABS(gExtItr.max(binPhi), phiMax, 1e-6);

  // Create a second Extent
  Extent gExtCopy;
  gExtCopy.extend(gExt);

  CHECK_CLOSE_ABS(gExtCopy.min(binX), 15_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.max(binX), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.min(binY), -3_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.max(binY), 3_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.min(binZ), -10_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.max(binZ), 10_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.min(binR), rMin, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.max(binR), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.min(binPhi), phiMin, 1e-6);
  CHECK_CLOSE_ABS(gExtCopy.max(binPhi), phiMax, 1e-6);

  // Check containment
  Extent unbound;
  BOOST_CHECK(unbound.contains(gExt));
  BOOST_CHECK(unbound.contains(gExtCopy));

  // Check application of an envelope on it
  ExtentEnvelope xEnvelopes = zeroEnvelopes;
  xEnvelopes[binX] = {1., 2.};

  // Take the extent and extend by an envelope
  Extent envelope(xEnvelopes);
  gExt.extend(envelope);
  // Changed ones
  CHECK_CLOSE_ABS(gExt.min(binX), 14_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binX), 20_mm, 1e-6);
  // Unchanged ones
  CHECK_CLOSE_ABS(gExt.min(binY), -3_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binY), 3_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.min(binZ), -10_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binZ), 10_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.min(binR), rMin, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binR), 18_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.min(binPhi), phiMin, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binPhi), phiMax, 1e-6);

  // Fill it with envelope
  Extent gExtEnv(envelope);
  gExtEnv.extend(vertices.begin(), vertices.end());
  // Changed ones
  CHECK_CLOSE_ABS(gExtEnv.min(binX), 14_mm, 1e-6);
  CHECK_CLOSE_ABS(gExtEnv.max(binX), 20_mm, 1e-6);

  // Check the set method
  gExt.set(binX, 2_mm, 8_mm);
  CHECK_CLOSE_ABS(gExt.min(binX), 2_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binX), 8_mm, 1e-6);

  // Radius can not go below 0
  gExt.set(binR, -2_mm, 18_mm);
  CHECK_CLOSE_ABS(gExt.min(binR), 0_mm, 1e-6);
  CHECK_CLOSE_ABS(gExt.max(binR), 18_mm, 1e-6);

  // Take an Extent and add a constraint
  Extent gExtConst;
  gExtConst.set(binR, 0., 5.);
  Extent gExtNonConst;
  BOOST_CHECK(!gExtNonConst.constrains(binR));
  gExtNonConst.addConstrain(gExtConst);
  BOOST_CHECK(gExtNonConst.constrains(binR));

  std::string tString = gExtConst.toString();
  BOOST_CHECK(!tString.empty());

  // Check single vertex containment
  Extent gExtVertexCheck;
  gExtVertexCheck.set(binR, 0., 5.);
  BOOST_CHECK(gExtVertexCheck.contains(Vector3(1., 0., 0.)));
  BOOST_CHECK(!gExtVertexCheck.contains(Vector3(6., 0., 0.)));
}

// Test that the constrains() check advances when the extend() method
// is used with a new binning type
BOOST_AUTO_TEST_CASE(ProtoSupportCaseTests) {
  std::vector<Vector3> vertices = {
      Vector3(15_mm, -3_mm, -10_mm), Vector3(18_mm, 0_mm, -10_mm),
      Vector3(15_mm, 3_mm, -10_mm),  Vector3(15_mm, -3_mm, 10_mm),
      Vector3(18_mm, 0_mm, 10_mm),   Vector3(15_mm, 3_mm, 10_mm)};

  Extent volumeExtent;
  volumeExtent.set(binZ, -300_mm, 300_mm);

  BOOST_CHECK(volumeExtent.constrains(binZ));
  BOOST_CHECK(!volumeExtent.constrains(binR));

  for (const auto& v : vertices) {
    volumeExtent.extend(v, {binR});
  }

  BOOST_CHECK(volumeExtent.constrains(binR));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
