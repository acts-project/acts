// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/PointBounds.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <stdexcept>
#include <vector>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

/// Unit test for creating compliant/non-compliant PointBounds object
BOOST_AUTO_TEST_CASE(PointBoundsConstruction) {
  double maxDistance = 3.;
  PointBounds pointBounds(maxDistance);
  BOOST_CHECK_EQUAL(pointBounds.type(), SurfaceBounds::ePoint);

  // Copy constructor
  PointBounds copied(pointBounds);
  BOOST_CHECK_EQUAL(copied.type(), SurfaceBounds::ePoint);
  BOOST_CHECK_EQUAL(copied.get(PointBounds::eR), maxDistance);
}

/// Unit test for PointBounds properties
BOOST_AUTO_TEST_CASE(PointBoundsProperties) {
  double maxDistance = 3.;
  PointBounds pointBounds(maxDistance);

  BOOST_CHECK_EQUAL(pointBounds.get(PointBounds::eR), maxDistance);
  BOOST_CHECK(pointBounds.isCartesian());

  // values array
  std::vector<double> values = pointBounds.values();
  BOOST_CHECK_EQUAL(values.size(), 1u);
  BOOST_CHECK_EQUAL(values[0], maxDistance);

  // center is the origin
  BOOST_CHECK_EQUAL(pointBounds.center(), Vector2::Zero());
}

/// Unit test for the inside() disc semantics (Cartesian x^2 + y^2 <= R^2)
BOOST_AUTO_TEST_CASE(PointBoundsInside) {
  double maxDistance = 3.;
  PointBounds pointBounds(maxDistance);

  BOOST_CHECK(pointBounds.inside(Vector2(0., 0.)));
  BOOST_CHECK(pointBounds.inside(Vector2(2., 2.)));      // r ~ 2.83 < 3
  BOOST_CHECK(pointBounds.inside(Vector2(3., 0.)));      // on boundary
  BOOST_CHECK(pointBounds.inside(Vector2(0., -3.)));     // on boundary
  BOOST_CHECK(!pointBounds.inside(Vector2(3., 3.)));     // r ~ 4.24 > 3
  BOOST_CHECK(!pointBounds.inside(Vector2(-2.5, -2.)));  // r ~ 3.2 > 3
}

/// Unit test for closestPoint (radial projection onto the boundary circle)
BOOST_AUTO_TEST_CASE(PointBoundsClosestPoint) {
  double maxDistance = 3.;
  PointBounds pointBounds(maxDistance);

  Vector2 closest =
      pointBounds.closestPoint(Vector2(6., 0.), SquareMatrix2::Identity());
  BOOST_CHECK_CLOSE(closest.x(), 3., 1e-9);
  BOOST_CHECK_SMALL(closest.y(), 1e-9);

  // at the origin any boundary point is valid; norm must equal R
  Vector2 fromOrigin =
      pointBounds.closestPoint(Vector2(0., 0.), SquareMatrix2::Identity());
  BOOST_CHECK_CLOSE(fromOrigin.norm(), 3., 1e-9);
}

/// Negative radius must throw
BOOST_AUTO_TEST_CASE(PointBoundsException) {
  BOOST_CHECK_THROW(PointBounds(-1.), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
