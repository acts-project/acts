// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE DiscTrapezoidal Bounds Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/DiscTrapezoidalBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit tests for DiscTrapezoidalBounds constrcuctors
BOOST_AUTO_TEST_CASE(DiscTrapezoidalBoundsConstruction) {
  double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0), averagePhi(0.0),
      stereo(0.1);
  // test default construction
  // DiscTrapezoidalBounds defaultConstructedDiscTrapezoidalBounds;  should be
  // deleted
  //
  /// Test construction with dimensions and default stereo
  BOOST_CHECK_EQUAL(
      DiscTrapezoidalBounds(minHalfX, maxHalfX, rMin, rMax, averagePhi).type(),
      SurfaceBounds::DiscTrapezoidal);
  //
  /// Test construction with all dimensions
  BOOST_CHECK_EQUAL(
      DiscTrapezoidalBounds(minHalfX, maxHalfX, rMin, rMax, averagePhi, stereo)
          .type(),
      SurfaceBounds::DiscTrapezoidal);
  //
  /// Copy constructor
  DiscTrapezoidalBounds original(minHalfX, maxHalfX, rMin, rMax, averagePhi);
  DiscTrapezoidalBounds copied(original);
  BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::DiscTrapezoidal);
}

/// Unit tests for DiscTrapezoidalBounds properties
BOOST_AUTO_TEST_CASE(DiscTrapezoidalBoundsProperties) {
  double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0),
      averagePhi(0.0) /*, stereo(0.1)*/;
  /// Test clone
  DiscTrapezoidalBounds discTrapezoidalBoundsObject(minHalfX, maxHalfX, rMin,
                                                    rMax, averagePhi);
  auto pClonedDiscTrapezoidalBounds = discTrapezoidalBoundsObject.clone();
  BOOST_CHECK_NE(pClonedDiscTrapezoidalBounds, nullptr);
  delete pClonedDiscTrapezoidalBounds;
  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(discTrapezoidalBoundsObject.type(),
                    SurfaceBounds::DiscTrapezoidal);
  //
  /// Test distanceToBoundary
  Vector2D origin(0., 0.);
  Vector2D outside(30., 0.);
  Vector2D inSurface(2., 0.0);
  CHECK_CLOSE_REL(discTrapezoidalBoundsObject.distanceToBoundary(origin), 2.0,
                  1e-6);
  CHECK_CLOSE_REL(discTrapezoidalBoundsObject.distanceToBoundary(outside), 24.0,
                  1e-6);
  //
  /// Test dump
  boost::test_tools::output_test_stream dumpOuput;
  discTrapezoidalBoundsObject.toStream(dumpOuput);
  BOOST_CHECK(dumpOuput.is_equal(
      "Acts::DiscTrapezoidalBounds:  (innerRadius, outerRadius, hMinX, "
      "hMaxX, hlengthY, hPhiSector, averagePhi, rCenter, stereo) = "
      "(2.0000000, 6.0000000, 1.0000000, 5.0000000, 0.7922870, 0.9851108, "
      "0.0000000, 2.5243378, 0.0000000)"));
  //
  /// Test inside
  BOOST_CHECK(
      discTrapezoidalBoundsObject.inside(inSurface, BoundaryCheck(true)));
  BOOST_CHECK(
      !discTrapezoidalBoundsObject.inside(outside, BoundaryCheck(true)));
  //
  /// Test rMin
  CHECK_CLOSE_REL(discTrapezoidalBoundsObject.rMin(), rMin, 1e-6);
  //
  /// Test rMax
  CHECK_CLOSE_REL(discTrapezoidalBoundsObject.rMax(), rMax, 1e-6);
  //
  /// Test averagePhi
  CHECK_SMALL(discTrapezoidalBoundsObject.averagePhi(), 1e-9);
  //
  /// Test rCenter (redundant; not configurable)
  CHECK_CLOSE_REL(discTrapezoidalBoundsObject.rCenter(), 2.524337798, 1e-6);
  //
  /// Test halfPhiSector (redundant; not configurable)
  CHECK_SMALL(discTrapezoidalBoundsObject.stereo(), 1e-6);
  //
  /// Test minHalflengthX
  CHECK_CLOSE_REL(discTrapezoidalBoundsObject.minHalflengthX(), minHalfX, 1e-6);
  //
  /// Test maxHalflengthX
  CHECK_CLOSE_REL(discTrapezoidalBoundsObject.maxHalflengthX(), maxHalfX, 1e-6);
  //
  /// Test halflengthY
  CHECK_CLOSE_REL(discTrapezoidalBoundsObject.halflengthY(), 0.792286991, 1e-6);
}
/// Unit test for testing DiscTrapezoidalBounds assignment
BOOST_AUTO_TEST_CASE(DiscTrapezoidalBoundsAssignment) {
  double minHalfX(1.0), maxHalfX(5.0), rMin(2.0), rMax(6.0), averagePhi(0.0),
      stereo(0.1);
  DiscTrapezoidalBounds discTrapezoidalBoundsObject(minHalfX, maxHalfX, rMin,
                                                    rMax, averagePhi, stereo);
  // operator == not implemented in this class
  //
  /// Test assignment
  DiscTrapezoidalBounds assignedDiscTrapezoidalBoundsObject(2.1, 6.6, 3.4, 4.2,
                                                            33.);
  assignedDiscTrapezoidalBoundsObject = discTrapezoidalBoundsObject;
  BOOST_CHECK_EQUAL(assignedDiscTrapezoidalBoundsObject,
                    discTrapezoidalBoundsObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
