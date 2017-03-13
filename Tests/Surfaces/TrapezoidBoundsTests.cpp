// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Trapezoid Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank

#include <boost/test/data/test_case.hpp>
// leave blank

#include <boost/test/output_test_stream.hpp>
// leave blank

//
#include "ACTS/Surfaces/TrapezoidBounds.hpp"
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
  /// Unit test for creating compliant/non-compliant TrapezoidBounds object
  BOOST_AUTO_TEST_CASE(TrapezoidBoundsConstruction)
  {
    double minHalfX(1.), maxHalfX(6.), halfY(2.);
    double alpha(std::atan(5. / 4.)), beta(std::atan(1.));
    //
    // default construction  deleted
    // TrapezoidBounds defaultConstructedTrapezoidBounds;
    //
    /// Test construction with defining half lengths
    BOOST_TEST(TrapezoidBounds(minHalfX, maxHalfX, halfY).type()
               == SurfaceBounds::Trapezoid);
    //
    /// Test construction with lengths and angles
    BOOST_TEST(TrapezoidBounds(minHalfX, halfY, alpha, beta).type()
               == SurfaceBounds::Trapezoid);
    /// Copy constructor
    TrapezoidBounds original(minHalfX, maxHalfX, halfY);
    TrapezoidBounds copied(original);
    BOOST_TEST(copied.type() == SurfaceBounds::Trapezoid);
  }

  /// Unit tests for TrapezoidBounds properties
  BOOST_AUTO_TEST_CASE(TrapezoidBoundsProperties, *utf::expected_failures(3))
  {
    double minHalfX(1.), maxHalfX(6.), halfY(2.);
    double alpha(std::atan(5. / 4.)), beta(std::atan(1.));
    //
    TrapezoidBounds trapezoidBoundsObject(minHalfX, maxHalfX, halfY);
    TrapezoidBounds asymmetricTrapezoidBoundsObject(
        minHalfX, halfY, alpha, beta);
    /// Test clone
    auto pClonedTrapezoidBounds = trapezoidBoundsObject.clone();
    BOOST_TEST(bool(pClonedTrapezoidBounds));
    delete pClonedTrapezoidBounds;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_TEST(trapezoidBoundsObject.type() == SurfaceBounds::Trapezoid);
    //
    /// Test minHalflengthX
    BOOST_TEST(trapezoidBoundsObject.minHalflengthX() == minHalfX);
    //
    /// Test maxHalfLengthX
    BOOST_TEST(trapezoidBoundsObject.maxHalflengthX() == maxHalfX);
    //
    /// Test maxHalfLengthX for asymmetric trapezoid, fails
    BOOST_TEST(asymmetricTrapezoidBoundsObject.maxHalflengthX() == maxHalfX);
    //
    /// Test halflengthY
    BOOST_TEST(trapezoidBoundsObject.halflengthY() == halfY);
    //
    /// Test alpha for symmetric trapezoid, will fail
    BOOST_TEST(trapezoidBoundsObject.alpha() == std::atan(5. / 4.));
    //
    /// Test alpha for asymmetric trapezoid
    BOOST_TEST(asymmetricTrapezoidBoundsObject.alpha() == alpha);
    //
    /// Test beta for symmetric trapezoid, will fail
    BOOST_TEST(trapezoidBoundsObject.beta() == std::atan(5. / 4.));
    //
    /// Test beta for asymmetric trapezoid
    BOOST_TEST(asymmetricTrapezoidBoundsObject.beta() == beta);
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outside(30., 0.);
    Vector2D inRectangle(2., 0.5);
    BOOST_TEST(trapezoidBoundsObject.distanceToBoundary(origin)
               == -2.);  // makes sense
    BOOST_TEST(trapezoidBoundsObject.distanceToBoundary(outside)
               == std::hypot(2., 24.));  // ok
    //
    /// Test vertices
    std::vector<Vector2D> expectedVertices{
        {1., -2.}, {6., 2.}, {-6., 2.}, {-1., -2.}};
    BOOST_TEST(trapezoidBoundsObject.vertices() == expectedVertices);
    /**
    for (auto i: trapezoidBoundsObject.vertices()){
      std::cout<<i[0]<<", "<<i[1]<<std::endl;
    }**/
    //
    /// Test boundingBox
    BOOST_TEST(trapezoidBoundsObject.boundingBox() == RectangleBounds(6., 2.));
    //

    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    trapezoidBoundsObject.dump(dumpOuput);
    BOOST_TEST(
        dumpOuput.is_equal("Acts::TrapezoidBounds:  (minHlenghtX, maxHlengthX, "
                           "hlengthY) = (1.0000000, 6.0000000, 2.0000000)"));
    //
    /// Test inside
    BOOST_TEST(trapezoidBoundsObject.inside(inRectangle, BoundaryCheck(true))
               == true);
    BOOST_TEST(trapezoidBoundsObject.inside(outside, BoundaryCheck(true))
               == false);
    //
    /// Test insideLoc0
    BOOST_TEST(trapezoidBoundsObject.insideLoc0(inRectangle) == true);
    //
    /// Test insideLoc1
    BOOST_TEST(trapezoidBoundsObject.insideLoc1(inRectangle) == true);
    //
  }
  /// Unit test for testing TrapezoidBounds assignment
  BOOST_AUTO_TEST_CASE(TrapezoidBoundsAssignment)
  {
    double          minHalfX(1.), maxHalfX(6.), halfY(2.);
    TrapezoidBounds trapezoidBoundsObject(minHalfX, maxHalfX, halfY);
    // operator == not implemented in this class
    //
    /// Test assignment
    TrapezoidBounds assignedTrapezoidBoundsObject(
        NaN, NaN, NaN, NaN);  // invalid object, in some sense
    assignedTrapezoidBoundsObject = trapezoidBoundsObject;
    BOOST_TEST(assignedTrapezoidBoundsObject == trapezoidBoundsObject);
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts
