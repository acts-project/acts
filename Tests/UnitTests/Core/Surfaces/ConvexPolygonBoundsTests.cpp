// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <chrono>
#include <iostream>
#include <memory>

#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

using vec2 = Acts::Vector2D;
template <int N>
using poly = Acts::ConvexPolygonBounds<N>;

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Surfaces)

BOOST_AUTO_TEST_CASE(convexity_test) {
  std::vector<vec2> vertices;
  vertices = {{0, 0}, {1, 0}, {0.2, 0.2}, {0, 1}};
  { BOOST_CHECK_THROW(poly<4> quad(vertices), AssertionFailureException); }

  vertices = {{0, 0}, {1, 0}, {0.8, 0.8}, {0, 1}};
  {
    // wrong number of vertices
    BOOST_CHECK_THROW(poly<3> trip{vertices}, AssertionFailureException);
  }
  {
    poly<4> quad = {vertices};
    BOOST_CHECK(quad.convex());
  }

  // this one is self intersecting
  vertices = {{0, 0}, {1, 0}, {0.5, 1}, {0.9, 1.2}};
  { BOOST_CHECK_THROW(poly<4> quad{vertices}, AssertionFailureException); }

  // this one is not
  vertices = {{0, 0}, {1, 0}, {0.9, 1.2}, {0.5, 1}};
  {
    poly<4> quad = {vertices};
    BOOST_CHECK(quad.convex());
  }

  vertices = {{0, 0}, {1, 0}, {0.8, 0.5}, {1, 1}, {0, 1}};
  { BOOST_CHECK_THROW(poly<5> pent(vertices), AssertionFailureException); }

  vertices = {{0, 0}, {1, 0}, {1.1, 0.5}, {1, 1}, {0, 1}};
  {
    poly<5> pent{vertices};
    BOOST_CHECK(pent.convex());
  }
}

BOOST_AUTO_TEST_CASE(construction_test) {
  std::vector<vec2> vertices;

  // triangle
  vertices = {{0, 0}, {1, 0}, {0.5, 1}};
  poly<3> triangle(vertices);

  RectangleBounds bb = triangle.boundingBox();
  BOOST_CHECK_EQUAL(bb.min(), Vector2D(0, 0));
  BOOST_CHECK_EQUAL(bb.max(), Vector2D(1., 1));

  BoundaryCheck bc(true);

  BOOST_CHECK(triangle.inside({0.2, 0.2}, bc));
  BOOST_CHECK(!triangle.inside({0.4, 0.9}, bc));
  BOOST_CHECK(!triangle.inside({0.8, 0.8}, bc));
  BOOST_CHECK(!triangle.inside({0.3, -0.2}, bc));

  CHECK_CLOSE_ABS(triangle.distanceToBoundary({0.2, 0.2}), -0.0894427, 1e-6);
  CHECK_CLOSE_ABS(triangle.distanceToBoundary({0.4, 0.9}), 0.0447213, 1e-6);
  CHECK_CLOSE_ABS(triangle.distanceToBoundary({0.8, 0.8}), 0.1788854, 1e-6);
  CHECK_CLOSE_ABS(triangle.distanceToBoundary({0.3, -0.2}), 0.2, 1e-6);

  // rectangular poly
  vertices = {{0, 0}, {1, 0}, {0.9, 1.2}, {0.5, 1}};
  poly<4> quad(vertices);

  bb = quad.boundingBox();
  BOOST_CHECK_EQUAL(bb.min(), Vector2D(0, 0));
  BOOST_CHECK_EQUAL(bb.max(), Vector2D(1, 1.2));

  BOOST_CHECK(quad.inside({0.2, 0.2}, bc));
  BOOST_CHECK(!quad.inside({0.4, 0.9}, bc));
  BOOST_CHECK(quad.inside({0.8, 0.8}, bc));
  BOOST_CHECK(!quad.inside({0.3, -0.2}, bc));

  CHECK_CLOSE_ABS(quad.distanceToBoundary({0.2, 0.2}), -0.089442, 1e-6);
  CHECK_CLOSE_ABS(quad.distanceToBoundary({0.4, 0.9}), 0.044721, 1e-6);
  CHECK_CLOSE_ABS(quad.distanceToBoundary({0.8, 0.8}), -0.132872, 1e-6);
  CHECK_CLOSE_ABS(quad.distanceToBoundary({0.3, -0.2}), 0.2, 1e-6);
}

BOOST_AUTO_TEST_CASE(construction_test_dynamic) {
  using poly = ConvexPolygonBounds<PolygonDynamic>;

  std::vector<vec2> vertices;

  // triangle
  vertices = {{0, 0}, {1, 0}, {0.5, 1}};
  poly triangle(vertices);

  RectangleBounds bb = triangle.boundingBox();
  BOOST_CHECK_EQUAL(bb.min(), Vector2D(0, 0));
  BOOST_CHECK_EQUAL(bb.max(), Vector2D(1., 1));

  BoundaryCheck bc(true);

  BOOST_CHECK(triangle.inside({0.2, 0.2}, bc));
  BOOST_CHECK(!triangle.inside({0.4, 0.9}, bc));
  BOOST_CHECK(!triangle.inside({0.8, 0.8}, bc));
  BOOST_CHECK(!triangle.inside({0.3, -0.2}, bc));

  CHECK_CLOSE_ABS(triangle.distanceToBoundary({0.2, 0.2}), -0.0894427, 1e-6);
  CHECK_CLOSE_ABS(triangle.distanceToBoundary({0.4, 0.9}), 0.0447213, 1e-6);
  CHECK_CLOSE_ABS(triangle.distanceToBoundary({0.8, 0.8}), 0.1788854, 1e-6);
  CHECK_CLOSE_ABS(triangle.distanceToBoundary({0.3, -0.2}), 0.2, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace Acts
