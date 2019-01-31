// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Trapezoid Bounds Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant TrapezoidBounds object
  BOOST_AUTO_TEST_CASE(TrapezoidBoundsConstruction)
  {
    double minHalfX(1.), maxHalfX(6.), halfY(2.);
    //
    // default construction  deleted
    // TrapezoidBounds defaultConstructedTrapezoidBounds;
    //
    /// Test construction with defining half lengths
    BOOST_CHECK_EQUAL(TrapezoidBounds(minHalfX, maxHalfX, halfY).type(),
                      SurfaceBounds::Trapezoid);
    /// Copy constructor
    TrapezoidBounds original(minHalfX, maxHalfX, halfY);
    TrapezoidBounds copied(original);
    BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::Trapezoid);
  }

  /// Unit tests for TrapezoidBounds properties
  BOOST_AUTO_TEST_CASE(TrapezoidBoundsProperties, *utf::expected_failures(3))
  {
    double minHalfX(1.), maxHalfX(6.), halfY(2.);
    //
    TrapezoidBounds trapezoidBoundsObject(minHalfX, maxHalfX, halfY);
    /// Test clone
    auto pClonedTrapezoidBounds = trapezoidBoundsObject.clone();
    BOOST_CHECK_NE(pClonedTrapezoidBounds, nullptr);
    delete pClonedTrapezoidBounds;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_CHECK_EQUAL(trapezoidBoundsObject.type(), SurfaceBounds::Trapezoid);
    //
    /// Test minHalflengthX
    BOOST_CHECK_EQUAL(trapezoidBoundsObject.minHalflengthX(), minHalfX);
    //
    /// Test maxHalfLengthX
    BOOST_CHECK_EQUAL(trapezoidBoundsObject.maxHalflengthX(), maxHalfX);
    //
    /// Test halflengthY
    BOOST_CHECK_EQUAL(trapezoidBoundsObject.halflengthY(), halfY);
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outside(30., 0.);
    Vector2D inRectangle(2., 0.5);
    CHECK_CLOSE_REL(trapezoidBoundsObject.distanceToBoundary(origin),
                    -2.,
                    1e-6);  // makes sense
    CHECK_CLOSE_REL(trapezoidBoundsObject.distanceToBoundary(outside),
                    std::hypot(2., 24.),
                    1e-6);  // ok
    //
    /// Test vertices
    std::vector<Vector2D> expectedVertices{
        {1., -2.}, {6., 2.}, {-6., 2.}, {-1., -2.}};
    const auto& actualVertices = trapezoidBoundsObject.vertices();
    BOOST_CHECK_EQUAL_COLLECTIONS(actualVertices.cbegin(),
                                  actualVertices.cend(),
                                  expectedVertices.cbegin(),
                                  expectedVertices.cend());
    /**
    for (auto i: trapezoidBoundsObject.vertices()){
      std::cout<<i[0]<<", "<<i[1]<<std::endl;
    }**/
    //
    /// Test boundingBox
    BOOST_CHECK_EQUAL(trapezoidBoundsObject.boundingBox(),
                      RectangleBounds(6., 2.));
    //

    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    trapezoidBoundsObject.dump(dumpOuput);
    BOOST_CHECK(
        dumpOuput.is_equal("Acts::TrapezoidBounds:  (minHlengthX, maxHlengthX, "
                           "hlengthY) = (1.0000000, 6.0000000, 2.0000000)"));
    //
    /// Test inside
    BOOST_CHECK(trapezoidBoundsObject.inside(inRectangle, BoundaryCheck(true)));
    BOOST_CHECK(!trapezoidBoundsObject.inside(outside, BoundaryCheck(true)));
  }
  /// Unit test for testing TrapezoidBounds assignment
  BOOST_AUTO_TEST_CASE(TrapezoidBoundsAssignment)
  {
    double          minHalfX(1.), maxHalfX(6.), halfY(2.);
    TrapezoidBounds trapezoidBoundsObject(minHalfX, maxHalfX, halfY);
    // operator == not implemented in this class
    //
    /// Test assignment
    TrapezoidBounds assignedTrapezoidBoundsObject(10., 20., 14.2);
    assignedTrapezoidBoundsObject = trapezoidBoundsObject;
    BOOST_CHECK_EQUAL(assignedTrapezoidBoundsObject, trapezoidBoundsObject);
  }

  BOOST_AUTO_TEST_CASE(TrapezoidBounds_toVariantData)
  {
    double          minHlX = 10;
    double          maxHlX = 15;
    double          hlY    = 5;
    TrapezoidBounds trap(minHlX, maxHlX, hlY);
    variant_data    var_data = trap.toVariantData();

    std::cout << var_data << std::endl;

    variant_map var_map = boost::get<variant_map>(var_data);
    BOOST_CHECK_EQUAL(var_map.get<std::string>("type"), "TrapezoidBounds");
    variant_map pl = var_map.get<variant_map>("payload");
    BOOST_CHECK_EQUAL(pl.get<double>("minHalfX"), minHlX);
    BOOST_CHECK_EQUAL(pl.get<double>("maxHalfX"), maxHlX);
    BOOST_CHECK_EQUAL(pl.get<double>("halfY"), hlY);

    TrapezoidBounds trap2(var_data);
    BOOST_CHECK_EQUAL(trap.minHalflengthX(), trap2.minHalflengthX());
    BOOST_CHECK_EQUAL(trap.maxHalflengthX(), trap2.maxHalflengthX());
    BOOST_CHECK_EQUAL(trap.halflengthY(), trap2.halflengthY());
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
