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

#include "Acts/Surfaces/EllipseBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit test for creating compliant/non-compliant EllipseBounds object
BOOST_AUTO_TEST_CASE(EllipseBoundsConstruction) {
  double minRad0(10.), maxRad0(15.), minRad1(25.), maxRad1(30.),
      phiSector(M_PI / 2.), averagePhi(0.);
  // test default construction
  // EllipseBounds defaultConstructedEllipseBounds;  //deleted
  //
  /// Test construction with dimensions
  BOOST_CHECK_EQUAL(
      EllipseBounds(minRad0, maxRad0, minRad1, maxRad1, phiSector, averagePhi)
          .type(),
      SurfaceBounds::eEllipse);
  //
  /// Copy constructor
  EllipseBounds original(minRad0, maxRad0, minRad1, maxRad1, phiSector,
                         averagePhi);
  EllipseBounds copied(original);
  BOOST_CHECK_EQUAL(copied, original);
}

// Streaning and recreation test
BOOST_AUTO_TEST_CASE(EllipseBoundsRecreation) {
  double minRad0(10.), maxRad0(15.), minRad1(25.), maxRad1(30.),
      phiSector(M_PI / 2.), averagePhi(0.);
  // const bool symmetric(false);
  EllipseBounds original(minRad0, maxRad0, minRad1, maxRad1, phiSector,
                         averagePhi);
  auto valvector = original.values();
  std::array<double, EllipseBounds::eSize> values;
  std::copy_n(valvector.begin(), EllipseBounds::eSize, values.begin());
  EllipseBounds recreated(values);
  BOOST_CHECK_EQUAL(recreated, original);
}

// Unit tests for AnnulusBounds exception throwing
BOOST_AUTO_TEST_CASE(ConeBoundsExceptions) {
  double minRad0(10.), maxRad0(15.), minRad1(25.), maxRad1(30.),
      phiSector(M_PI / 2.), averagePhi(0.);
  // Exception for opening minR0 < 0
  BOOST_CHECK_THROW(
      EllipseBounds(-minRad0, maxRad0, minRad1, maxRad1, phiSector, averagePhi),
      std::logic_error);
  // Exception for opening maxR0 < 0
  BOOST_CHECK_THROW(
      EllipseBounds(minRad0, -maxRad0, minRad1, maxRad1, phiSector, averagePhi),
      std::logic_error);
  // Exception for opening minR0 and maxR0 < 0
  BOOST_CHECK_THROW(EllipseBounds(-minRad0, -maxRad0, minRad1, maxRad1,
                                  phiSector, averagePhi),
                    std::logic_error);
  // Exception for swapped minR0/maxR0
  BOOST_CHECK_THROW(
      EllipseBounds(maxRad0, minRad0, minRad1, maxRad1, phiSector, averagePhi),
      std::logic_error);
  // Exception for opening minR1 < 0
  BOOST_CHECK_THROW(
      EllipseBounds(minRad0, maxRad0, -minRad1, maxRad1, phiSector, averagePhi),
      std::logic_error);
  // Exception for opening maxR1 < 0
  BOOST_CHECK_THROW(
      EllipseBounds(minRad0, maxRad0, minRad1, -maxRad1, phiSector, averagePhi),
      std::logic_error);
  // Exception for opening maxR1 < 0
  BOOST_CHECK_THROW(EllipseBounds(minRad0, maxRad0, -minRad1, -maxRad1,
                                  phiSector, averagePhi),
                    std::logic_error);
  // Exception for swapped minR1/maxR1
  BOOST_CHECK_THROW(
      EllipseBounds(minRad0, maxRad0, maxRad1, minRad1, phiSector, averagePhi),
      std::logic_error);
  // Exception for negative phiSector
  BOOST_CHECK_THROW(
      EllipseBounds(minRad0, maxRad0, minRad1, maxRad1, -phiSector, averagePhi),
      std::logic_error);
  // Exception for average phi out of bound
  BOOST_CHECK_THROW(
      EllipseBounds(minRad0, maxRad0, minRad1, maxRad1, phiSector, 4.),
      std::logic_error);
}

/// Unit tests for EllipseBounds properties
BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(CylinderBoundsProperties, 1)
BOOST_AUTO_TEST_CASE(EllipseBoundsProperties) {
  double minRad0(10.), minRad1(15.), maxRad0(15.), maxRad1(20.), averagePhi(0.),
      phiSector(M_PI / 2.);
  /// Test clone
  EllipseBounds ellipseBoundsObject(minRad0, maxRad0, minRad1, maxRad1,
                                    phiSector, averagePhi);
  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(ellipseBoundsObject.type(), SurfaceBounds::eEllipse);
  //
  // clone already tested
  //
  /// Test distanceToBoundary
  Vector2D origin(0., 0.);
  Vector2D outsideBy15(0., 30.);
  Vector2D inRectangle(17., 11.);
  CHECK_CLOSE_REL(ellipseBoundsObject.distanceToBoundary(origin), 10.,
                  1e-6);  // makes sense
  CHECK_CLOSE_REL(ellipseBoundsObject.distanceToBoundary(outsideBy15), 15.,
                  1e-6);  // fails, not clear why
  //
  /// Test rMinX
  BOOST_CHECK_EQUAL(ellipseBoundsObject.get(EllipseBounds::eMinR0), minRad0);
  //
  /// Test rMinY
  BOOST_CHECK_EQUAL(ellipseBoundsObject.get(EllipseBounds::eMinR1), minRad1);
  //
  /// Test rMaxX
  BOOST_CHECK_EQUAL(ellipseBoundsObject.get(EllipseBounds::eMaxR0), maxRad0);
  //
  /// Test rMaxY
  BOOST_CHECK_EQUAL(ellipseBoundsObject.get(EllipseBounds::eMaxR1), maxRad1);
  //
  /// Test averagePhi
  BOOST_CHECK_EQUAL(ellipseBoundsObject.get(EllipseBounds::eAveragePhi),
                    averagePhi);
  //
  /// Test vertices
  // std::vector<Vector2D> expectedVertices{{15, 0}, {0, 20}, {-15, 0}, {0,
  // -20}}; const auto& actualVertices = ellipseBoundsObject.vertices(4);
  // BOOST_CHECK_EQUAL_COLLECTIONS(actualVertices.cbegin(),
  // actualVertices.cend(),
  //                              expectedVertices.cbegin(),
  //                              expectedVertices.cend());
  //
  /// Test boundingBox
  BOOST_CHECK_EQUAL(ellipseBoundsObject.boundingBox(),
                    RectangleBounds(15., 20.));
  //
  /// Test halfPhiSector
  BOOST_CHECK_EQUAL(ellipseBoundsObject.get(EllipseBounds::eHalfPhiSector),
                    M_PI / 2.);
  //
  /// Test dump
  boost::test_tools::output_test_stream dumpOuput;
  ellipseBoundsObject.toStream(dumpOuput);
  BOOST_CHECK(dumpOuput.is_equal(
      "Acts::EllipseBounds:  (innerRadius0, outerRadius0, innerRadius1, "
      "outerRadius1, hPhiSector, averagePhi) = (10.0000000, 15.0000000, "
      "15.0000000, "
      "20.0000000, 0.0000000, 1.5707963, 0.0000000)"));
  //
  /// Test inside
  BOOST_CHECK(!ellipseBoundsObject.inside(inRectangle, BoundaryCheck(true)));
  // dont understand why this is so:
  BOOST_CHECK(!ellipseBoundsObject.inside(outsideBy15, BoundaryCheck(true)));
}
/// Unit test for testing EllipseBounds assignment
BOOST_AUTO_TEST_CASE(EllipseBoundsAssignment) {
  double minRad0(10.), minRad1(15.), maxRad0(15.), maxRad1(20.), averagePhi(0.),
      phiSector(M_PI / 2.);
  EllipseBounds ellipseBoundsObject(minRad0, minRad1, maxRad0, maxRad1,
                                    averagePhi, phiSector);
  EllipseBounds similarlyConstructeEllipseBoundsObject(
      minRad0, minRad1, maxRad0, maxRad1, averagePhi, phiSector);
  /// Test operator ==
  BOOST_CHECK_EQUAL(ellipseBoundsObject,
                    similarlyConstructeEllipseBoundsObject);
  //
  /// Test assignment
  EllipseBounds assignedEllipseBoundsObject(11., 12., 17., 18., 1.);
  // object, in some sense
  assignedEllipseBoundsObject = ellipseBoundsObject;
  BOOST_CHECK_EQUAL(assignedEllipseBoundsObject, ellipseBoundsObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
