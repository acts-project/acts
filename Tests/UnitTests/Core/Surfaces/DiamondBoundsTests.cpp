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

#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit test for creating compliant/non-compliant DiamondBounds object
BOOST_AUTO_TEST_CASE(DiamondBoundsConstruction) {
  double minHalfX(10.), midHalfX(20.), maxHalfX(15.), halfY1(5.), halfY2(7.);
  // test default construction
  // DiamondBounds defaultConstructedDiamondBounds;  //deleted
  //
  /// Test construction with dimensions
  // DiamondBounds d(minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
  BOOST_CHECK_EQUAL(
      DiamondBounds(minHalfX, midHalfX, maxHalfX, halfY1, halfY2).type(),
      SurfaceBounds::Diamond);
  //
  /// Copy constructor
  DiamondBounds original(minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
  DiamondBounds copied(original);
  BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::Diamond);

  // invalid inputs
  BOOST_CHECK_THROW(
      DiamondBounds db(midHalfX, minHalfX, maxHalfX, halfY1, halfY2),
      AssertionFailureException);
  BOOST_CHECK_THROW(
      DiamondBounds db(minHalfX, maxHalfX, midHalfX, halfY1, halfY2),
      AssertionFailureException);
}
/// Unit tests for DiamondBounds properties
BOOST_AUTO_TEST_CASE(DiamondBoundsProperties) {
  double minHalfX(10.), midHalfX(50.), maxHalfX(30.), halfY1(10.), halfY2(20.);
  /// Test clone
  DiamondBounds diamondBoundsObject(minHalfX, midHalfX, maxHalfX, halfY1,
                                    halfY2);
  auto pClonedDiamondBoundsObject = diamondBoundsObject.clone();
  BOOST_CHECK_NE(pClonedDiamondBoundsObject, nullptr);
  delete pClonedDiamondBoundsObject;
  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(diamondBoundsObject.type(), SurfaceBounds::Diamond);
  // //redundant test
  //
  /// Test minHalflengthX() NOTE: Naming violation
  BOOST_CHECK_EQUAL(diamondBoundsObject.x1(), minHalfX);
  //
  /// Test medHalflengthX() NOTE: Naming violation
  BOOST_CHECK_EQUAL(diamondBoundsObject.x2(), midHalfX);
  //
  /// Test maxHalflengthX() NOTE: Naming violation
  BOOST_CHECK_EQUAL(diamondBoundsObject.x3(), maxHalfX);
  //
  /// Test halflengthY1() NOTE: Naming violation
  BOOST_CHECK_EQUAL(diamondBoundsObject.y1(), halfY1);
  //
  /// Test halflengthY2() NOTE: Naming violation
  BOOST_CHECK_EQUAL(diamondBoundsObject.y2(), halfY2);
  //
  /// Test boundingBox
  BOOST_CHECK_EQUAL(diamondBoundsObject.boundingBox(),
                    RectangleBounds(50., 20.));
  //
  // clone already tested
  //
  /// Test distanceToBoundary
  Vector2D origin(0., 0.);
  Vector2D outsideBy10(0., 30.);
  Vector2D inRectangle(15., 0.);
  CHECK_CLOSE_REL(diamondBoundsObject.distanceToBoundary(origin), -10.,
                  1e-6);  // makes sense
  CHECK_CLOSE_REL(diamondBoundsObject.distanceToBoundary(outsideBy10), 10.,
                  1e-6);  // ok
  //
  /// Test dump
  // Acts::DiamondBounds:  (minHlengthX, medHlengthX, maxHlengthX, hlengthY1,
  // hlengthY2 ) = (30.0000000, 10.0000000, 50.0000000, 10.0000000,
  // 20.0000000)
  diamondBoundsObject.toStream(std::cout);
  boost::test_tools::output_test_stream dumpOuput;
  diamondBoundsObject.toStream(dumpOuput);
  BOOST_CHECK(
      dumpOuput.is_equal("Acts::DiamondBounds:  (x1, x2, "
                         "x3, y1, y2 ) = (10.0000000, "
                         "50.0000000, 30.0000000, 10.0000000, 20.0000000)"));
  //
  /// Test inside
  BOOST_CHECK(diamondBoundsObject.inside(origin, BoundaryCheck(true)));
  // dont understand why this is so:
  BOOST_CHECK(!diamondBoundsObject.inside(outsideBy10, BoundaryCheck(true)));
  //
  /// Test vertices (does this need to be implemented in this class??
  // auto v=diamondBoundsObject.vertices();
  std::vector<Vector2D> referenceVertices{
      {-minHalfX, -halfY1}, {minHalfX, -halfY1}, {midHalfX, 0.},
      {maxHalfX, halfY2},   {-maxHalfX, halfY2}, {-midHalfX, 0.}};
  const auto& actualVertices = diamondBoundsObject.vertices();
  BOOST_CHECK_EQUAL_COLLECTIONS(actualVertices.cbegin(), actualVertices.cend(),
                                referenceVertices.cbegin(),
                                referenceVertices.cend());
}
/// Unit test for testing DiamondBounds assignment
BOOST_AUTO_TEST_CASE(DiamondBoundsAssignment) {
  double minHalfX(10.), midHalfX(20.), maxHalfX(15.), halfY1(5.), halfY2(7.);
  DiamondBounds diamondBoundsObject(minHalfX, midHalfX, maxHalfX, halfY1,
                                    halfY2);
  DiamondBounds similarlyConstructeDiamondBoundsObject(
      minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
  /// Test operator ==
  BOOST_CHECK_EQUAL(diamondBoundsObject,
                    similarlyConstructeDiamondBoundsObject);
  //
  /// Test assignment
  DiamondBounds assignedDiamondBoundsObject(0, 0, 0, 0, 0);  // invalid
  // object, in some sense
  assignedDiamondBoundsObject = diamondBoundsObject;
  BOOST_CHECK_EQUAL(assignedDiamondBoundsObject, diamondBoundsObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
