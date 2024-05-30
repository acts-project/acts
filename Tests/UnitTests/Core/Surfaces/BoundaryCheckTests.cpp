// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "BoundaryCheckTestsRefs.hpp"

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Surfaces)
// See: https://en.wikipedia.org/wiki/Bounding_volume
//
// Aligned box w/ simple check
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxSimple) {
  BoundaryCheck check(true);
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  BOOST_CHECK(check.isInside({0, 0}, ll, ur));
  BOOST_CHECK(!check.isInside({2, 2}, ll, ur));
  BOOST_CHECK(!check.isInside({0, 2}, ll, ur));
  BOOST_CHECK(!check.isInside({2, 0}, ll, ur));
}
// Aligned box w/ tolerance check along first axis
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxToleranceLoc0) {
  BoundaryCheck check(true, false, 1.5, 0.0);
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  BOOST_CHECK(check.isInside({0, 0}, ll, ur));
  BOOST_CHECK(check.isInside({2, 2}, ll, ur));
  BOOST_CHECK(!check.isInside({4, 4}, ll, ur));
  BOOST_CHECK(check.isInside({0, 2}, ll, ur));
  BOOST_CHECK(check.isInside({2, 0}, ll, ur));
}

BOOST_AUTO_TEST_CASE(BoundaryCheckBoxDistance) {
  BoundaryCheck bcheck(true);

  for (std::size_t i = 0; i < rectTestPoints.size(); i++) {
    const Vector2& testPoint = rectTestPoints.at(i);
    double refDistance = rectDistances.at(i);
    Vector2 ll(rectDimensions.xmin, rectDimensions.ymin);
    Vector2 ur(rectDimensions.xmax, rectDimensions.ymax);
    double distance = bcheck.distance(testPoint, ll, ur);
    CHECK_CLOSE_REL(refDistance, distance, 1e-6);
  }

  for (std::size_t i = 0; i < rectShiftedTestPoints.size(); i++) {
    const Vector2& testPoint = rectShiftedTestPoints.at(i);
    double refDistance = rectShiftedDistances.at(i);
    Vector2 ll(rectShiftedDimensions.xmin, rectShiftedDimensions.ymin);
    Vector2 ur(rectShiftedDimensions.xmax, rectShiftedDimensions.ymax);
    double distance = bcheck.distance(testPoint, ll, ur);
    CHECK_CLOSE_REL(refDistance, distance, 1e-6);
  }
}

// Aligned box w/ covariance check
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxCovariance) {
  SquareMatrix2 cov;
  cov << 1, 0.5, 0.5, 2;
  BoundaryCheck check(cov, 3.0);
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  BOOST_CHECK(check.isInside({0, 0}, ll, ur));
  BOOST_CHECK(check.isInside({2, 2}, ll, ur));
  BOOST_CHECK(!check.isInside({4, 4}, ll, ur));
  BOOST_CHECK(check.isInside({0, 3}, ll, ur));
  BOOST_CHECK(check.isInside({3, 0}, ll, ur));
}

BOOST_AUTO_TEST_CASE(BoundaryCheckPolyDistance) {
  // we check a box again, but this time described as a poly

  BoundaryCheck bcheck(true);

  for (std::size_t i = 0; i < rectTestPoints.size(); i++) {
    const Vector2& testPoint = rectTestPoints.at(i);
    double refDistance = rectDistances.at(i);
    double distance = bcheck.distance(testPoint, rectVertices);
    CHECK_CLOSE_REL(refDistance, distance, 1e-6);
  }

  for (std::size_t i = 0; i < rectShiftedTestPoints.size(); i++) {
    const Vector2& testPoint = rectShiftedTestPoints.at(i);
    double refDistance = rectShiftedDistances.at(i);
    double distance = bcheck.distance(testPoint, rectShiftedVertices);
    CHECK_CLOSE_REL(refDistance, distance, 1e-6);
  }
}

// Triangle w/ simple check
BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleSimple) {
  Vector2 vertices[] = {{-2, 0}, {2, 0}, {0, 2}};
  BoundaryCheck check(true);
  BOOST_CHECK(check.isInside({0, 0}, vertices));
  BOOST_CHECK(check.isInside({0, 1}, vertices));
  BOOST_CHECK(!check.isInside({2, 2}, vertices));
  BOOST_CHECK(!check.isInside({0, -1}, vertices));
}
// Triangle w/ covariance check
BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleCovariance) {
  Vector2 vertices[] = {{-2, 0}, {2, 0}, {0, 2}};
  SquareMatrix2 cov;
  cov << 0.5, 0, 0, 0.5;
  BoundaryCheck check(cov, 4.1);
  BOOST_CHECK(check.isInside({0, 0}, vertices));
  BOOST_CHECK(check.isInside({0, 1}, vertices));
  BOOST_CHECK(check.isInside({0, 2}, vertices));
  BOOST_CHECK(check.isInside({0, 3}, vertices));
  BOOST_CHECK(check.isInside({0, 4}, vertices));
  BOOST_CHECK(!check.isInside({0, 5}, vertices));
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
