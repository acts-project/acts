// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/detail/IntersectionHelper2D.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <limits>

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Surfaces)

void basicChecks(bool circleCase = false) {
  double rY = 10.;
  double rX = circleCase ? rY : 5.;

  // Line along x - not intersecting
  Vector2D start(12., 0.);
  Vector2D direction(0., -1.);

  auto nosol = circleCase ? detail::IntersectionHelper2D::intersectCircle(
                                rY, start, direction)
                          : detail::IntersectionHelper2D::intersectEllipse(
                                rX, rY, start, direction);
  BOOST_CHECK(not nosol.first);
  BOOST_CHECK(not nosol.second);

  start = Vector2D(4., -4.);
  auto twosol = circleCase ? detail::IntersectionHelper2D::intersectCircle(
                                 rY, start, direction)
                           : detail::IntersectionHelper2D::intersectEllipse(
                                 rX, rY, start, direction);

  BOOST_CHECK(twosol.first);
  BOOST_CHECK(twosol.second);

  start = Vector2D(-4., 10.);
  direction = Vector2D(1., 0.);

  auto onesolY = circleCase ? detail::IntersectionHelper2D::intersectCircle(
                                  rY, start, direction)
                            : detail::IntersectionHelper2D::intersectEllipse(
                                  rX, rY, start, direction);

  BOOST_CHECK(onesolY.first);
  CHECK_CLOSE_ABS(onesolY.first.position.x(), 0., s_epsilon);
  BOOST_CHECK(not onesolY.second);

  start = Vector2D(rX, -4);
  direction = Vector2D(0., 1.);

  auto onesolX = circleCase ? detail::IntersectionHelper2D::intersectCircle(
                                  rY, start, direction)
                            : detail::IntersectionHelper2D::intersectEllipse(
                                  rX, rY, start, direction);

  BOOST_CHECK(onesolX.first);
  CHECK_CLOSE_ABS(onesolX.first.position.y(), 0., s_epsilon);
  BOOST_CHECK(not onesolX.second);
}

/// Unit test for creating Ellipse intersection
BOOST_AUTO_TEST_CASE(LineLineIntersection) {
  Vector2D start(1., 1.);
  Vector2D dir(1., 1.);

  // Not possible
  auto solution = detail::IntersectionHelper2D::intersectSegment(
      Vector2D(5., 3.), Vector2D(6., 4), start, dir);

  BOOST_CHECK(not solution);

  // Possible
  solution = detail::IntersectionHelper2D::intersectSegment(
      Vector2D(5., 3.), Vector2D(3., -1.), start, dir);

  BOOST_CHECK(solution);
}

/// Unit test for creating Ellipse intersection
BOOST_AUTO_TEST_CASE(EllipseIntersection) {
  // Basic checks for ellipse
  basicChecks();

  // Specific checks for ellipse
  double radiusX = 450.;
  double radiusY = 275.;

  Vector2D start(-500., -300.);
  Vector2D direction = Vector2D(10., 4.).normalized();

  auto solution = detail::IntersectionHelper2D::intersectEllipse(
      radiusX, radiusY, start, direction);

  // Numerically checked / per hand calculated
  BOOST_CHECK(solution.first);

  CHECK_CLOSE_ABS(solution.first.position.x(), -283.68, 0.01);
  CHECK_CLOSE_ABS(solution.first.position.y(), -213.47, 0.01);
  BOOST_CHECK(solution.first.pathLength > 0.);

  BOOST_CHECK(solution.second);

  CHECK_CLOSE_ABS(solution.second.position.x(), 433.65, 0.01);
  CHECK_CLOSE_ABS(solution.second.position.y(), 73.46, 0.01);
  BOOST_CHECK(solution.second.pathLength > 0.);

  // Reverse checks will be done with circle (same code)
}

/// Unit test for creating Circle intersection
BOOST_AUTO_TEST_CASE(CircleIntersection) {
  // Basic checks for circle
  basicChecks(true);

  // Specific checks for circle
  double radius = 275.;

  Vector2D start(-500., -300.);
  Vector2D direction = Vector2D(1., 1.).normalized();

  auto solution =
      detail::IntersectionHelper2D::intersectCircle(radius, start, direction);

  // Numerically checked / per hand calculated
  BOOST_CHECK(solution.first);

  CHECK_CLOSE_ABS(solution.first.position.x(), -266.771, 0.001);
  CHECK_CLOSE_ABS(solution.first.position.y(), -66.771, 0.001);
  BOOST_CHECK(solution.first.pathLength > 0.);

  BOOST_CHECK(solution.second);

  CHECK_CLOSE_ABS(solution.second.position.x(), 66.771, 0.001);
  CHECK_CLOSE_ABS(solution.second.position.y(), 266.771, 0.001);
  BOOST_CHECK(solution.second.pathLength > 0.);

  // Reverse
  start = Vector2D(1500., 1700.);
  direction = Vector2D(1., 1.).normalized();
  solution =
      detail::IntersectionHelper2D::intersectCircle(radius, start, direction);

  BOOST_CHECK(solution.first);
  CHECK_CLOSE_ABS(solution.first.position.x(), 66.771, 0.001);
  CHECK_CLOSE_ABS(solution.first.position.y(), 266.771, 0.001);
  BOOST_CHECK(solution.first.pathLength < 0.);

  BOOST_CHECK(solution.second);
  CHECK_CLOSE_ABS(solution.second.position.x(), -266.771, 0.001);
  CHECK_CLOSE_ABS(solution.second.position.y(), -66.771, 0.001);
  BOOST_CHECK(solution.second.pathLength < 0.);

  // Reverse with reverse direction
  direction = Vector2D(-1., -1.).normalized();
  solution =
      detail::IntersectionHelper2D::intersectCircle(radius, start, direction);

  BOOST_CHECK(solution.first);
  CHECK_CLOSE_ABS(solution.first.position.x(), 66.771, 0.001);
  CHECK_CLOSE_ABS(solution.first.position.y(), 266.771, 0.001);
  BOOST_CHECK(solution.first.pathLength > 0.);

  BOOST_CHECK(solution.second);
  CHECK_CLOSE_ABS(solution.second.position.x(), -266.771, 0.001);
  CHECK_CLOSE_ABS(solution.second.position.y(), -66.771, 0.001);
  BOOST_CHECK(solution.second.pathLength > 0.);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts