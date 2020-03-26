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

#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit test for creating compliant/non-compliant TrapezoidBounds object
BOOST_AUTO_TEST_CASE(TrapezoidBoundsConstruction) {
  double minHalfX(1.), maxHalfX(6.), halfY(2.);
  //
  // default construction  deleted
  // TrapezoidBounds defaultConstructedTrapezoidBounds;
  //
  /// Test construction with defining half lengths
  BOOST_CHECK_EQUAL(TrapezoidBounds(minHalfX, maxHalfX, halfY).type(),
                    SurfaceBounds::eTrapezoid);
  /// Copy constructor
  TrapezoidBounds original(minHalfX, maxHalfX, halfY);
  TrapezoidBounds copied(original);
  BOOST_CHECK_EQUAL(copied, original);
}

/// Unit test for creating compliant/non-compliant TrapezoidBounds object
BOOST_AUTO_TEST_CASE(TrapezoidBoundsRecreated) {
  double minHalfX(1.), maxHalfX(6.), halfY(2.);
  /// Copy constructor
  TrapezoidBounds original(minHalfX, maxHalfX, halfY);
  // const bool symmetric(false);
  auto valvector = original.values();
  std::array<double, TrapezoidBounds::eSize> values;
  std::copy_n(valvector.begin(), TrapezoidBounds::eSize, values.begin());
  TrapezoidBounds recreated(values);
  BOOST_CHECK_EQUAL(original, recreated);
}

// Exception tests
BOOST_AUTO_TEST_CASE(TrapezoidBoundsException) {
  double minHalfX(1.), maxHalfX(6.), halfY(2.);

  // Negative x at min y
  BOOST_CHECK_THROW(TrapezoidBounds(-minHalfX, maxHalfX, halfY),
                    std::logic_error);

  // Negative x at max y
  BOOST_CHECK_THROW(TrapezoidBounds(minHalfX, -maxHalfX, halfY),
                    std::logic_error);

  // Negative x at miny and max y
  BOOST_CHECK_THROW(TrapezoidBounds(-minHalfX, -maxHalfX, halfY),
                    std::logic_error);

  // Negative y
  BOOST_CHECK_THROW(TrapezoidBounds(minHalfX, maxHalfX, -halfY),
                    std::logic_error);
}

/// Unit tests for TrapezoidBounds properties
BOOST_AUTO_TEST_CASE(TrapezoidBoundsProperties, *utf::expected_failures(3)) {
  double minHalfX(1.), maxHalfX(6.), halfY(2.);
  //
  TrapezoidBounds trapezoidBoundsObject(minHalfX, maxHalfX, halfY);
  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(trapezoidBoundsObject.type(), SurfaceBounds::eTrapezoid);
  //
  /// Test minHalflengthX
  BOOST_CHECK_EQUAL(
      trapezoidBoundsObject.get(TrapezoidBounds::eHalfLengthXnegY), minHalfX);
  //
  /// Test maxHalfLengthX
  BOOST_CHECK_EQUAL(
      trapezoidBoundsObject.get(TrapezoidBounds::eHalfLengthXposY), maxHalfX);
  //
  /// Test halflengthY
  BOOST_CHECK_EQUAL(trapezoidBoundsObject.get(TrapezoidBounds::eHalfLengthY),
                    halfY);
  //
  /// Test distanceToBoundary
  Vector2D origin(0., 0.);
  Vector2D outside(30., 0.);
  Vector2D inRectangle(2., 0.5);
  CHECK_CLOSE_REL(trapezoidBoundsObject.distanceToBoundary(origin), -2.,
                  1e-6);  // makes sense
  CHECK_CLOSE_REL(trapezoidBoundsObject.distanceToBoundary(outside),
                  std::hypot(2., 24.),
                  1e-6);  // ok
  //
  /// Test vertices
  std::vector<Vector2D> expectedVertices{
      {1., -2.}, {6., 2.}, {-6., 2.}, {-1., -2.}};
  const auto& actualVertices = trapezoidBoundsObject.vertices();
  BOOST_CHECK_EQUAL_COLLECTIONS(actualVertices.cbegin(), actualVertices.cend(),
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
  trapezoidBoundsObject.toStream(dumpOuput);
  BOOST_CHECK(
      dumpOuput.is_equal("Acts::TrapezoidBounds:  (minHlengthX, maxHlengthX, "
                         "hlengthY) = (1.0000000, 6.0000000, 2.0000000)"));
  //
  /// Test inside
  BOOST_CHECK(trapezoidBoundsObject.inside(inRectangle, BoundaryCheck(true)));
  BOOST_CHECK(!trapezoidBoundsObject.inside(outside, BoundaryCheck(true)));
}
/// Unit test for testing TrapezoidBounds assignment
BOOST_AUTO_TEST_CASE(TrapezoidBoundsAssignment) {
  double minHalfX(1.), maxHalfX(6.), halfY(2.);
  TrapezoidBounds trapezoidBoundsObject(minHalfX, maxHalfX, halfY);
  // operator == not implemented in this class
  //
  /// Test assignment
  TrapezoidBounds assignedTrapezoidBoundsObject(10., 20., 14.2);
  assignedTrapezoidBoundsObject = trapezoidBoundsObject;
  BOOST_CHECK_EQUAL(assignedTrapezoidBoundsObject, trapezoidBoundsObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
