// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Diamond Bounds Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/DiamondBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant DiamondBounds object
  BOOST_AUTO_TEST_CASE(DiamondBoundsConstruction)
  {
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
  BOOST_AUTO_TEST_CASE(DiamondBoundsProperties)
  {
    double minHalfX(10.), midHalfX(50.), maxHalfX(30.), halfY1(10.),
        halfY2(20.);
    /// Test clone
    DiamondBounds diamondBoundsObject(
        minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
    auto pClonedDiamondBoundsObject = diamondBoundsObject.clone();
    BOOST_CHECK_NE(pClonedDiamondBoundsObject, nullptr);
    delete pClonedDiamondBoundsObject;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_CHECK_EQUAL(diamondBoundsObject.type(), SurfaceBounds::Diamond);
    // //redundant test
    //
    /// Test minHalflengthX() NOTE: Naming violation
    BOOST_CHECK_EQUAL(diamondBoundsObject.minHalflengthX(), minHalfX);
    //
    /// Test medHalflengthX() NOTE: Naming violation
    BOOST_CHECK_EQUAL(diamondBoundsObject.medHalflengthX(), midHalfX);
    //
    /// Test maxHalflengthX() NOTE: Naming violation
    BOOST_CHECK_EQUAL(diamondBoundsObject.maxHalflengthX(), maxHalfX);
    //
    /// Test halflengthY1() NOTE: Naming violation
    BOOST_CHECK_EQUAL(diamondBoundsObject.halflengthY1(), halfY1);
    //
    /// Test halflengthY2() NOTE: Naming violation
    BOOST_CHECK_EQUAL(diamondBoundsObject.halflengthY2(), halfY2);
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
    CHECK_CLOSE_REL(diamondBoundsObject.distanceToBoundary(origin),
                    -10.,
                    1e-6);  // makes sense
    CHECK_CLOSE_REL(diamondBoundsObject.distanceToBoundary(outsideBy10),
                    10.,
                    1e-6);  // ok
    //
    /// Test dump
    // Acts::DiamondBounds:  (minHlengthX, medHlengthX, maxHlengthX, hlengthY1,
    // hlengthY2 ) = (30.0000000, 10.0000000, 50.0000000, 10.0000000,
    // 20.0000000)
    diamondBoundsObject.dump(std::cout);
    boost::test_tools::output_test_stream dumpOuput;
    diamondBoundsObject.dump(dumpOuput);
    BOOST_CHECK(
        dumpOuput.is_equal("Acts::DiamondBounds:  (minHlengthX, medHlengthX, "
                           "maxHlengthX, hlengthY1, hlengthY2 ) = (10.0000000, "
                           "50.0000000, 30.0000000, 10.0000000, 20.0000000)"));
    //
    /// Test inside
    BOOST_CHECK(diamondBoundsObject.inside(origin, BoundaryCheck(true)));
    // dont understand why this is so:
    BOOST_CHECK(!diamondBoundsObject.inside(outsideBy10, BoundaryCheck(true)));
    //
    /// Test vertices (does this need to be implemented in this class??
    // auto v=diamondBoundsObject.vertices();
    std::vector<Vector2D> referenceVertices{{minHalfX, -halfY1},
                                            {midHalfX, 0.},
                                            {maxHalfX, halfY2},
                                            {-maxHalfX, halfY2},
                                            {-midHalfX, 0.},
                                            {-minHalfX, -halfY1}};
    const auto& actualVertices = diamondBoundsObject.vertices();
    BOOST_CHECK_EQUAL_COLLECTIONS(actualVertices.cbegin(),
                                  actualVertices.cend(),
                                  referenceVertices.cbegin(),
                                  referenceVertices.cend());
  }
  /// Unit test for testing DiamondBounds assignment
  BOOST_AUTO_TEST_CASE(DiamondBoundsAssignment)
  {
    double minHalfX(10.), midHalfX(20.), maxHalfX(15.), halfY1(5.), halfY2(7.);
    DiamondBounds diamondBoundsObject(
        minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
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

  BOOST_AUTO_TEST_CASE(DiamondBounds_toVariantData)
  {
    double minHalfX(10.), midHalfX(50.), maxHalfX(30.), halfY1(10.),
        halfY2(20.);
    DiamondBounds diam(minHalfX, midHalfX, maxHalfX, halfY1, halfY2);
    variant_data  var_data = diam.toVariantData();

    std::cout << var_data << std::endl;

    variant_map var_map = boost::get<variant_map>(var_data);
    BOOST_CHECK_EQUAL(var_map.get<std::string>("type"), "DiamondBounds");
    variant_map pl = var_map.get<variant_map>("payload");
    BOOST_CHECK_EQUAL(pl.get<double>("minHalfX"), minHalfX);
    BOOST_CHECK_EQUAL(pl.get<double>("medHalfX"), midHalfX);
    BOOST_CHECK_EQUAL(pl.get<double>("maxHalfX"), maxHalfX);
    BOOST_CHECK_EQUAL(pl.get<double>("minY"), halfY1);
    BOOST_CHECK_EQUAL(pl.get<double>("maxY"), halfY2);

    DiamondBounds diam2(var_data);
    BOOST_CHECK_EQUAL(diam.minHalflengthX(), diam2.minHalflengthX());
    BOOST_CHECK_EQUAL(diam.medHalflengthX(), diam2.medHalflengthX());
    BOOST_CHECK_EQUAL(diam.maxHalflengthX(), diam2.maxHalflengthX());
    BOOST_CHECK_EQUAL(diam.halflengthY1(), diam2.halflengthY1());
    BOOST_CHECK_EQUAL(diam.halflengthY2(), diam2.halflengthY2());
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
