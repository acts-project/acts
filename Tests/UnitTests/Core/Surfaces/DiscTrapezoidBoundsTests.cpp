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

#include <limits>

#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit tests for DiscTrapezoidBounds constrcuctors
BOOST_AUTO_TEST_CASE(DiscTrapezoidBoundsConstruction) {
  double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0), averagePhi(0.0),
      stereo(0.1);
  // test default construction
  // DiscTrapezoidBounds defaultConstructedDiscTrapezoidBounds;  should be
  // deleted
  //
  /// Test construction with dimensions and default stereo
  BOOST_CHECK_EQUAL(
      DiscTrapezoidBounds(minHalfX, maxHalfX, rMin, rMax, averagePhi).type(),
      SurfaceBounds::DiscTrapezoidal);
  //
  /// Test construction with all dimensions
  BOOST_CHECK_EQUAL(
      DiscTrapezoidBounds(minHalfX, maxHalfX, rMin, rMax, averagePhi, stereo)
          .type(),
      SurfaceBounds::DiscTrapezoidal);
  //
  /// Copy constructor
  DiscTrapezoidBounds original(minHalfX, maxHalfX, rMin, rMax, averagePhi);
  DiscTrapezoidBounds copied(original);
  BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::DiscTrapezoidal);
}

/// Unit tests for DiscTrapezoidBounds properties
BOOST_AUTO_TEST_CASE(DiscTrapezoidBoundsProperties) {
  double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0),
      averagePhi(0.0) /*, stereo(0.1)*/;
  /// Test clone
  DiscTrapezoidBounds DiscTrapezoidBoundsObject(minHalfX, maxHalfX, rMin, rMax,
                                                averagePhi);
  auto pClonedDiscTrapezoidBounds = DiscTrapezoidBoundsObject.clone();
  BOOST_CHECK_NE(pClonedDiscTrapezoidBounds, nullptr);
  delete pClonedDiscTrapezoidBounds;
  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(DiscTrapezoidBoundsObject.type(),
                    SurfaceBounds::DiscTrapezoidal);
  //
  /// Test distanceToBoundary
  Vector2D origin(0., 0.);
  Vector2D outside(30., 0.);
  Vector2D inSurface(2., 0.0);
  CHECK_CLOSE_REL(DiscTrapezoidBoundsObject.distanceToBoundary(origin), 2.0,
                  1e-6);
  CHECK_CLOSE_REL(DiscTrapezoidBoundsObject.distanceToBoundary(outside), 24.0,
                  1e-6);
  //
  /// Test dump
  boost::test_tools::output_test_stream dumpOuput;
  DiscTrapezoidBoundsObject.toStream(dumpOuput);
  BOOST_CHECK(dumpOuput.is_equal(
      "Acts::DiscTrapezoidBounds:  (innerRadius, outerRadius, hMinX, "
      "hMaxX, hlengthY, hPhiSector, averagePhi, rCenter, stereo) = "
      "(2.0000000, 6.0000000, 1.0000000, 5.0000000, 0.7922870, 0.9851108, "
      "0.0000000, 2.5243378, 0.0000000)"));
  //
  /// Test inside
  BOOST_CHECK(DiscTrapezoidBoundsObject.inside(inSurface, BoundaryCheck(true)));
  BOOST_CHECK(!DiscTrapezoidBoundsObject.inside(outside, BoundaryCheck(true)));
  //
  /// Test rMin
  CHECK_CLOSE_REL(DiscTrapezoidBoundsObject.rMin(), rMin, 1e-6);
  //
  /// Test rMax
  CHECK_CLOSE_REL(DiscTrapezoidBoundsObject.rMax(), rMax, 1e-6);
  //
  /// Test averagePhi
  CHECK_SMALL(DiscTrapezoidBoundsObject.averagePhi(), 1e-9);
  //
  /// Test rCenter (redundant; not configurable)
  CHECK_CLOSE_REL(DiscTrapezoidBoundsObject.rCenter(), 2.524337798, 1e-6);
  //
  /// Test halfPhiSector (redundant; not configurable)
  CHECK_SMALL(DiscTrapezoidBoundsObject.stereo(), 1e-6);
  //
  /// Test minHalflengthX
  CHECK_CLOSE_REL(DiscTrapezoidBoundsObject.minHalflengthX(), minHalfX, 1e-6);
  //
  /// Test maxHalflengthX
  CHECK_CLOSE_REL(DiscTrapezoidBoundsObject.maxHalflengthX(), maxHalfX, 1e-6);
  //
  /// Test halflengthY
  CHECK_CLOSE_REL(DiscTrapezoidBoundsObject.halflengthY(), 0.792286991, 1e-6);
}
/// Unit test for testing DiscTrapezoidBounds assignment
BOOST_AUTO_TEST_CASE(DiscTrapezoidBoundsAssignment) {
  double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0), averagePhi(0.0),
      stereo(0.1);
  DiscTrapezoidBounds DiscTrapezoidBoundsObject(minHalfX, maxHalfX, rMin, rMax,
                                                averagePhi, stereo);
  // operator == not implemented in this class
  //
  /// Test assignment
  DiscTrapezoidBounds assignedDiscTrapezoidBoundsObject(2.1, 6.6, 3.4, 4.2,
                                                        33.);
  assignedDiscTrapezoidBoundsObject = DiscTrapezoidBoundsObject;
  BOOST_CHECK_EQUAL(assignedDiscTrapezoidBoundsObject,
                    DiscTrapezoidBoundsObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
