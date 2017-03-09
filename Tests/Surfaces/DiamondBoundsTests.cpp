// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Diamond Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank

#include <boost/test/data/test_case.hpp>
// leave blank

#include <boost/test/output_test_stream.hpp>
// leave blank

//
#include "ACTS/Surfaces/DiamondBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
//
#include <limits>

// namespace bdata = boost::unit_test::data;
namespace utf    = boost::unit_test;
const double inf = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces);
  /// Unit test for creating compliant/non-compliant DiamondBounds object
  BOOST_AUTO_TEST_CASE(DiamondBoundsConstruction)
  {
    double minHalfX(10.), midHalfX(15.), maxHalfX(20.), halfY1(5.), halfY2(7.);
    // test default construction
    // DiamondBounds defaultConstructedDiamondBounds;  //deleted
    //
    /// Test construction with dimensions
    // DiamondBounds d(minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
    BOOST_TEST(
        DiamondBounds(minHalfX, midHalfX, maxHalfX, halfY1, halfY2).type()
        == SurfaceBounds::Diamond);
    //
    /// Copy constructor
    DiamondBounds original(minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
    DiamondBounds copied(original);
    BOOST_TEST(copied.type() == SurfaceBounds::Diamond);
  }
  /// Unit tests for DiamondBounds properties
  BOOST_AUTO_TEST_CASE(DiamondBoundsProperties)
  {
    double minHalfX(30.), midHalfX(10.), maxHalfX(50.), halfY1(10.),
        halfY2(20.);
    /// Test clone
    DiamondBounds diamondBoundsObject(
        minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
    auto pClonedDiamondBoundsObject = diamondBoundsObject.clone();
    BOOST_TEST(bool(pClonedDiamondBoundsObject));
    delete pClonedDiamondBoundsObject;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_TEST(diamondBoundsObject.type() == SurfaceBounds::Diamond);
    // //redundant test
    //
    /// Test minHalflengthX() NOTE: Naming violation
    BOOST_TEST(diamondBoundsObject.minHalflengthX() == minHalfX);
    //
    /// Test medHalflengthX() NOTE: Naming violation
    BOOST_TEST(diamondBoundsObject.medHalflengthX() == midHalfX);
    //
    /// Test maxHalflengthX() NOTE: Naming violation
    BOOST_TEST(diamondBoundsObject.maxHalflengthX() == maxHalfX);
    //
    /// Test halflengthY1() NOTE: Naming violation
    BOOST_TEST(diamondBoundsObject.halflengthY1() == halfY1);
    //
    /// Test halflengthY2() NOTE: Naming violation
    BOOST_TEST(diamondBoundsObject.halflengthY2() == halfY2);
    //
    //
    BOOST_TEST_MESSAGE(
        "The following two tests pass but do not match the documentation:");
    /// Test alpha1()
    // double openingAngle1 = diamondBoundsObject.alpha1();
    BOOST_TEST(diamondBoundsObject.alpha1() == -M_PI / 4.);
    //
    /// Test alpha2()
    // double openingAngle2 = diamondBoundsObject.alpha2();
    BOOST_TEST(diamondBoundsObject.alpha2() == -M_PI / 4.);
    //
    /// Test boundingBox
    BOOST_TEST(diamondBoundsObject.boundingBox() == RectangleBounds(50., 20.));
    //
    // clone already tested
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outsideBy10(0., 30.);
    Vector2D inRectangle(15., 0.);
    BOOST_TEST(diamondBoundsObject.distanceToBoundary(origin)
               == -10.);  // makes sense
    BOOST_TEST(diamondBoundsObject.distanceToBoundary(outsideBy10)
               == 10.);  // ok
    //
    /// Test dump
    // Acts::DiamondBounds:  (minHlenghtX, medHlengthX, maxHlengthX, hlengthY1,
    // hlengthY2 ) = (30.0000000, 10.0000000, 50.0000000, 10.0000000,
    // 20.0000000)
    diamondBoundsObject.dump(std::cout);
    boost::test_tools::output_test_stream dumpOuput;
    diamondBoundsObject.dump(dumpOuput);
    BOOST_TEST(
        dumpOuput.is_equal("Acts::DiamondBounds:  (minHlengthX, medHlengthX, "
                           "maxHlengthX, hlengthY1, hlengthY2 ) = (30.0000000, "
                           "10.0000000, 50.0000000, 10.0000000, 20.0000000)"));
    //
    /// Test inside
    BOOST_TEST(diamondBoundsObject.inside(origin, BoundaryCheck(true)) == true);
    // dont understand why this is so:
    BOOST_TEST(diamondBoundsObject.inside(outsideBy10, BoundaryCheck(true))
               == true);
    //
    /// Test insideLoc0 (only checks for inside bounding rectangle...don't
    /// understand answer here)
    BOOST_TEST(diamondBoundsObject.insideLoc0(inRectangle) == false);
    //
    /// Test insideLoc1 (only checks for inside bounding rectangle)
    BOOST_TEST(diamondBoundsObject.insideLoc1(inRectangle) == true);
    //
    /// Test vertices (does this need to be implemented in this class??
    // auto v=diamondBoundsObject.vertices();
    std::vector<Vector2D> referenceVertices{{30., -10.},
                                            {10., 0.},
                                            {50., 20.},
                                            {-50., 20.},
                                            {-10., 0.},
                                            {-30., -10.}};
    BOOST_TEST(diamondBoundsObject.vertices() == referenceVertices);
  }
  /// Unit test for testing DiamondBounds assignment
  BOOST_AUTO_TEST_CASE(DiamondBoundsAssignment)
  {
    double minHalfX(10.), midHalfX(15.), maxHalfX(20.), halfY1(5.), halfY2(7.);
    DiamondBounds diamondBoundsObject(
        minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
    DiamondBounds similarlyConstructeDiamondBoundsObject(
        minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
    /// Test operator ==
    BOOST_TEST(diamondBoundsObject == similarlyConstructeDiamondBoundsObject);
    //
    /// Test assignment
    DiamondBounds assignedDiamondBoundsObject(
        NaN, NaN, NaN, NaN, NaN);  // invalid
    // object, in some sense
    assignedDiamondBoundsObject = diamondBoundsObject;
    BOOST_TEST(assignedDiamondBoundsObject == diamondBoundsObject);
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts
