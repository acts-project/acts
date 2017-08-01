// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Triangle Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank

#include <boost/test/data/test_case.hpp>
// leave blank

#include <boost/test/output_test_stream.hpp>
// leave blank

#include <limits>

#include "ACTS/Surfaces/TriangleBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace utf    = boost::unit_test;
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces);
  /// Unit test for creating compliant/non-compliant TriangleBounds object
  BOOST_AUTO_TEST_CASE(TriangleBoundsConstruction)
  {
    std::vector<Vector2D> vertices{
        {1., 1.}, {4., 1.}, {4., 5.}};  // 3-4-5 triangle
    // test default construction
    // TriangleBounds defaultConstructedTriangleBounds;  //deleted
    //
    /// Test construction with vertices
    BOOST_TEST(TriangleBounds(vertices).type() == SurfaceBounds::Triangle);
    //
    /// Copy constructor
    TriangleBounds original(vertices);
    TriangleBounds copied(original);
    BOOST_TEST(copied.type() == SurfaceBounds::Triangle);
  }

  /// Unit tests for TriangleBounds properties
  BOOST_AUTO_TEST_CASE(TriangleBoundsProperties, *utf::expected_failures(2))
  {
    std::vector<Vector2D> vertices{
        {1., 1.}, {4., 1.}, {4., 5.}};  // 3-4-5 triangle
    /// Test clone
    TriangleBounds triangleBoundsObject(vertices);
    auto           pClonedTriangleBounds = triangleBoundsObject.clone();
    BOOST_TEST(bool(pClonedTriangleBounds));
    delete pClonedTriangleBounds;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_TEST(triangleBoundsObject.type() == SurfaceBounds::Triangle);
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outside(30., 1.);
    Vector2D inTriangle(2., 1.5);
    BOOST_TEST(triangleBoundsObject.distanceToBoundary(origin)
               == std::sqrt(2.));  // makes sense
    BOOST_TEST(triangleBoundsObject.distanceToBoundary(outside) == 26.);  // ok
    //
    /// Test vertices : fail; there are in fact 6 vertices (10.03.2017)
    std::vector<Vector2D> expectedVertices(vertices);
    BOOST_TEST_MESSAGE(
        "Following two tests fail because the triangle has six vertices");
    BOOST_TEST(triangleBoundsObject.vertices().size() == (size_t)3);
    BOOST_TEST(triangleBoundsObject.vertices() == expectedVertices);
    // for (auto i: triangleBoundsObject.vertices()){
    //  std::cout<<i[0]<<", "<<i[1]<<std::endl;
    //}
    //
    /// Test boundingBox NOTE: Bounding box too big
    BOOST_TEST(triangleBoundsObject.boundingBox() == RectangleBounds(4., 5.));
    //

    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    triangleBoundsObject.dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal(
        "Acts::TriangleBounds:  generating vertices (X, Y)(1.0000000 , 1.0000000) \n\
(4.0000000 , 1.0000000) \n\
(4.0000000 , 5.0000000) "));
    //
    /// Test inside
    BOOST_TEST(triangleBoundsObject.inside(inTriangle, BoundaryCheck(true))
               == true);
    BOOST_TEST(triangleBoundsObject.inside(outside, BoundaryCheck(true))
               == false);
  }
  /// Unit test for testing TriangleBounds assignment
  BOOST_AUTO_TEST_CASE(TriangleBoundsAssignment)
  {
    std::vector<Vector2D> vertices{
        {1., 1.}, {4., 1.}, {4., 5.}};  // 3-4-5 triangle
    std::vector<Vector2D> invalid{{NaN, NaN}, {NaN, NaN}, {NaN, NaN}};
    TriangleBounds        triangleBoundsObject(vertices);
    // operator == not implemented in this class
    //
    /// Test assignment
    TriangleBounds assignedTriangleBoundsObject(
        invalid);  // invalid object, in some sense
    assignedTriangleBoundsObject = triangleBoundsObject;
    BOOST_TEST(assignedTriangleBoundsObject == triangleBoundsObject);
  }
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts
