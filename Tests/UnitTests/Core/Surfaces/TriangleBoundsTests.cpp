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

#include "Acts/Surfaces/TriangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit test for creating compliant/non-compliant TriangleBounds object
BOOST_AUTO_TEST_CASE(TriangleBoundsConstruction) {
  std::array<Vector2D, 3> vertices({{Vector2D(1., 1.), Vector2D(4., 1.),
                                     Vector2D(4., 5.)}});  // 3-4-5 triangle
  // test default construction
  // TriangleBounds defaultConstructedTriangleBounds;  //deleted
  //
  /// Test construction with vertices
  BOOST_CHECK_EQUAL(TriangleBounds(vertices).type(), SurfaceBounds::eTriangle);
  //
  /// Copy constructor
  TriangleBounds original(vertices);
  TriangleBounds copied(original);
  BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::eTriangle);
}

/// Unit tests for TriangleBounds properties
BOOST_AUTO_TEST_CASE(TriangleBoundsProperties) {
  std::array<Vector2D, 3> vertices({{Vector2D(1., 1.), Vector2D(4., 1.),
                                     Vector2D(4., 5.)}});  // 3-4-5 triangle
  /// Test clone
  TriangleBounds triangleBoundsObject(vertices);
  auto pClonedTriangleBounds = triangleBoundsObject.clone();
  BOOST_CHECK_NE(pClonedTriangleBounds, nullptr);
  delete pClonedTriangleBounds;
  //
  /// Test type() (redundant; already used in constructor confirmation)
  BOOST_CHECK_EQUAL(triangleBoundsObject.type(), SurfaceBounds::eTriangle);
  //
  /// Test distanceToBoundary
  Vector2D origin(0., 0.);
  Vector2D outside(30., 1.);
  Vector2D inTriangle(2., 1.5);
  CHECK_CLOSE_REL(triangleBoundsObject.distanceToBoundary(origin),
                  std::sqrt(2.),
                  1e-6);  // makes sense
  CHECK_CLOSE_REL(triangleBoundsObject.distanceToBoundary(outside), 26.,
                  1e-6);  // ok
  //
  /// Test vertices : fail; there are in fact 6 vertices (10.03.2017)
  std::array<Vector2D, 3> expectedVertices = vertices;
  BOOST_TEST_MESSAGE(
      "Following two tests fail because the triangle has six vertices");
  BOOST_CHECK_EQUAL(triangleBoundsObject.vertices().size(), (size_t)3);
  for (size_t i = 0; i < 3; i++) {
    CHECK_CLOSE_REL(triangleBoundsObject.vertices().at(i),
                    expectedVertices.at(i), 1e-6);
  }
  /// Test boundingBox NOTE: Bounding box too big
  BOOST_CHECK_EQUAL(triangleBoundsObject.boundingBox(),
                    RectangleBounds(4., 5.));
  //

  //
  /// Test dump
  boost::test_tools::output_test_stream dumpOuput;
  triangleBoundsObject.toStream(dumpOuput);
  BOOST_CHECK(dumpOuput.is_equal(
      "Acts::TriangleBounds:  generating vertices (X, Y)(1.0000000 , 1.0000000) \n\
(4.0000000 , 1.0000000) \n\
(4.0000000 , 5.0000000) "));
  //
  /// Test inside
  BOOST_CHECK(triangleBoundsObject.inside(inTriangle, BoundaryCheck(true)));
  BOOST_CHECK(!triangleBoundsObject.inside(outside, BoundaryCheck(true)));
}
/// Unit test for testing TriangleBounds assignment
BOOST_AUTO_TEST_CASE(TriangleBoundsAssignment) {
  std::array<Vector2D, 3> vertices({{Vector2D(1., 1.), Vector2D(4., 1.),
                                     Vector2D(4., 5)}});  // 3-4-5 triangle
  std::array<Vector2D, 3> invalid(
      {{Vector2D(-1, -1), Vector2D(-1, -1), Vector2D(-1, -1)}});
  TriangleBounds triangleBoundsObject(vertices);
  // operator == not implemented in this class
  //
  /// Test assignment
  TriangleBounds assignedTriangleBoundsObject(
      invalid);  // invalid object, in some sense
  assignedTriangleBoundsObject = triangleBoundsObject;
  const auto& assignedVertices = assignedTriangleBoundsObject.vertices();
  const auto& originalVertices = triangleBoundsObject.vertices();
  BOOST_CHECK_EQUAL_COLLECTIONS(
      assignedVertices.cbegin(), assignedVertices.cend(),
      originalVertices.cbegin(), originalVertices.cend());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
