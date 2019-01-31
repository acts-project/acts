// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Triangle Bounds Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/TriangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant TriangleBounds object
  BOOST_AUTO_TEST_CASE(TriangleBoundsConstruction)
  {
    std::array<Vector2D, 3> vertices({{Vector2D(1., 1.),
                                       Vector2D(4., 1.),
                                       Vector2D(4., 5.)}});  // 3-4-5 triangle
    // test default construction
    // TriangleBounds defaultConstructedTriangleBounds;  //deleted
    //
    /// Test construction with vertices
    BOOST_CHECK_EQUAL(TriangleBounds(vertices).type(), SurfaceBounds::Triangle);
    //
    /// Copy constructor
    TriangleBounds original(vertices);
    TriangleBounds copied(original);
    BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::Triangle);
  }

  /// Unit tests for TriangleBounds properties
  BOOST_AUTO_TEST_CASE(TriangleBoundsProperties)
  {
    std::array<Vector2D, 3> vertices({{Vector2D(1., 1.),
                                       Vector2D(4., 1.),
                                       Vector2D(4., 5.)}});  // 3-4-5 triangle
    /// Test clone
    TriangleBounds triangleBoundsObject(vertices);
    auto           pClonedTriangleBounds = triangleBoundsObject.clone();
    BOOST_CHECK_NE(pClonedTriangleBounds, nullptr);
    delete pClonedTriangleBounds;
    //
    /// Test type() (redundant; already used in constructor confirmation)
    BOOST_CHECK_EQUAL(triangleBoundsObject.type(), SurfaceBounds::Triangle);
    //
    /// Test distanceToBoundary
    Vector2D origin(0., 0.);
    Vector2D outside(30., 1.);
    Vector2D inTriangle(2., 1.5);
    CHECK_CLOSE_REL(triangleBoundsObject.distanceToBoundary(origin),
                    std::sqrt(2.),
                    1e-6);  // makes sense
    CHECK_CLOSE_REL(triangleBoundsObject.distanceToBoundary(outside),
                    26.,
                    1e-6);  // ok
    //
    /// Test vertices : fail; there are in fact 6 vertices (10.03.2017)
    std::array<Vector2D, 3> expectedVertices = vertices;
    BOOST_TEST_MESSAGE(
        "Following two tests fail because the triangle has six vertices");
    BOOST_CHECK_EQUAL(triangleBoundsObject.vertices().size(), (size_t)3);
    for (size_t i = 0; i < 3; i++) {
      CHECK_CLOSE_REL(
          triangleBoundsObject.vertices().at(i), expectedVertices.at(i), 1e-6);
    }
    /// Test boundingBox NOTE: Bounding box too big
    BOOST_CHECK_EQUAL(triangleBoundsObject.boundingBox(),
                      RectangleBounds(4., 5.));
    //

    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    triangleBoundsObject.dump(dumpOuput);
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
  BOOST_AUTO_TEST_CASE(TriangleBoundsAssignment)
  {
    std::array<Vector2D, 3> vertices({{Vector2D(1., 1.),
                                       Vector2D(4., 1.),
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
    BOOST_CHECK_EQUAL_COLLECTIONS(assignedVertices.cbegin(),
                                  assignedVertices.cend(),
                                  originalVertices.cbegin(),
                                  originalVertices.cend());
  }

  BOOST_AUTO_TEST_CASE(TriangleBounds_toVariantData)
  {
    std::array<Vector2D, 3> vertices({{Vector2D(1., 1.),
                                       Vector2D(4., 1.),
                                       Vector2D(4., 5.)}});  // 3-4-5 triangle
    TriangleBounds triangle(vertices);
    variant_data   var_data = triangle.toVariantData();

    std::cout << var_data << std::endl;

    variant_map var_map = boost::get<variant_map>(var_data);
    BOOST_CHECK_EQUAL(var_map.get<std::string>("type"), "TriangleBounds");
    variant_map pl = var_map.get<variant_map>("payload");

    variant_vector var_vertices = pl.get<variant_vector>("vertices");
    BOOST_CHECK_EQUAL(var_vertices.size(), 3);

    for (size_t i = 0; i < 3; i++) {
      Vector2D    exp = vertices.at(i);
      variant_map var = var_vertices.get<variant_map>(i);
      BOOST_CHECK_EQUAL(var.get<std::string>("type"), "Vector2D");
      variant_vector coords = var.get<variant_vector>("payload");

      BOOST_CHECK_EQUAL(exp.x(), coords.get<double>(0));
      BOOST_CHECK_EQUAL(exp.y(), coords.get<double>(1));
    }

    TriangleBounds triangle2(var_data);
    BOOST_CHECK_EQUAL(triangle2.vertices().size(), 3);
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
