// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/execution_monitor.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/ConvexPolygonBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Surfaces)

BOOST_AUTO_TEST_CASE(BoundaryToleranceConstructors) {
  using enum BoundaryTolerance::ToleranceMode;
  {
    // Test None constructor
    BoundaryTolerance tolerance = BoundaryTolerance::None();
    BOOST_CHECK(tolerance.toleranceMode() == None);
  }

  // Test AbsoluteEuclidean constructor
  {
    // Valid positive tolerance
    auto tolerance = BoundaryTolerance::AbsoluteEuclidean(1.0);
    BOOST_CHECK_EQUAL(tolerance.asAbsoluteEuclidean().tolerance, 1.0);
    BOOST_CHECK(tolerance.toleranceMode() == Extend);
    BOOST_CHECK(BoundaryTolerance::AbsoluteEuclidean(0.0).toleranceMode() ==
                None);

    // Valid negative tolerance
    tolerance = BoundaryTolerance::AbsoluteEuclidean(-1.0);
    BOOST_CHECK_EQUAL(tolerance.asAbsoluteEuclidean().tolerance, -1.0);
    BOOST_CHECK(tolerance.toleranceMode() == Shrink);
  }

  // Test Chi2Bound constructor
  {
    SquareMatrix2 cov;
    cov << 1, 0.5, 0.5, 2;

    // Valid positive chi2 bound
    auto tolerance = BoundaryTolerance::Chi2Bound(cov, 3.0);
    BOOST_CHECK_EQUAL(tolerance.asChi2Bound().maxChi2, 3.0);
    BOOST_CHECK(tolerance.toleranceMode() == Extend);
    BOOST_CHECK(BoundaryTolerance::Chi2Bound(cov, 0.0).toleranceMode() == None);

    // Valid negative chi2 bound
    tolerance = BoundaryTolerance::Chi2Bound(cov, -3.0);
    BOOST_CHECK_EQUAL(tolerance.asChi2Bound().maxChi2, -3.0);
    BOOST_CHECK(tolerance.toleranceMode() == Shrink);
  }

  // Test None constructor
  BoundaryTolerance::None();
}

// See: https://en.wikipedia.org/wiki/Bounding_volume
//
// Aligned box w/ simple check
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxSimple) {
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  RectangleBounds bounds(ll, ur);
  auto tolerance = BoundaryTolerance::None();
  BOOST_CHECK(bounds.inside({0, 0}, tolerance));
  BOOST_CHECK(!bounds.inside({2, 2}, tolerance));
  BOOST_CHECK(!bounds.inside({0, 2}, tolerance));
  BOOST_CHECK(!bounds.inside({2, 0}, tolerance));
}

// Aligned box w/ covariance check
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxCovariance) {
  SquareMatrix2 cov;
  cov << 1, 0.5, 0.5, 2;
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  RectangleBounds bounds(ll, ur);
  auto tolerance = BoundaryTolerance::Chi2Bound(cov.inverse(), 4);
  BOOST_CHECK(bounds.inside({0, 0}, tolerance));
  BOOST_CHECK(bounds.inside({2, 2}, tolerance));
  BOOST_CHECK(!bounds.inside({4, 4}, tolerance));
  BOOST_CHECK(bounds.inside({0, 3}, tolerance));
  BOOST_CHECK(bounds.inside({3, 0}, tolerance));
}

// Triangle w/ simple check
BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleSimple) {
  ConvexPolygonBounds<PolygonDynamic> poly({{-2, 0}, {2, 0}, {0, 2}});
  auto tolerance = BoundaryTolerance::None();
  BOOST_CHECK(poly.inside({0, 0}, tolerance));
  BOOST_CHECK(poly.inside({0, 1}, tolerance));
  BOOST_CHECK(!poly.inside({2, 2}, tolerance));
  BOOST_CHECK(!poly.inside({0, -1}, tolerance));
}
// Triangle w/ covariance check
BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleCovariance) {
  ConvexPolygonBounds<PolygonDynamic> poly({{-2, 0}, {2, 0}, {0, 2}});
  SquareMatrix2 cov;
  cov << 0.5, 0, 0, 0.5;
  auto tolerance = BoundaryTolerance::Chi2Bound(cov.inverse(), 8);
  BOOST_CHECK(poly.inside({0, 0}, tolerance));
  BOOST_CHECK(poly.inside({0, 1}, tolerance));
  BOOST_CHECK(poly.inside({0, 2}, tolerance));
  BOOST_CHECK(poly.inside({0, 3}, tolerance));
  BOOST_CHECK(poly.inside({0, 4}, tolerance));
  BOOST_CHECK(!poly.inside({0, 5}, tolerance));
}

BOOST_AUTO_TEST_CASE(BoundaryCheckDifferentTolerances) {
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  RectangleBounds bounds(ll, ur);

  {
    auto tolerance = BoundaryTolerance::None();
    BOOST_CHECK(bounds.inside({0, 0}, tolerance));
    BOOST_CHECK(!bounds.inside({2, 2}, tolerance));
    BOOST_CHECK(!bounds.inside({0, 2}, tolerance));
    BOOST_CHECK(!bounds.inside({2, 0}, tolerance));
  }

  {
    auto tolerance = BoundaryTolerance::Infinite();
    BOOST_CHECK(bounds.inside({0, 0}, tolerance));
    BOOST_CHECK(bounds.inside({2, 2}, tolerance));
    BOOST_CHECK(bounds.inside({0, 2}, tolerance));
    BOOST_CHECK(bounds.inside({2, 0}, tolerance));
  }

  {
    auto tolerance = BoundaryTolerance::AbsoluteEuclidean(1.1);
    BOOST_CHECK(bounds.inside({0, 0}, tolerance));
    BOOST_CHECK(!bounds.inside({2, 2}, tolerance));
    BOOST_CHECK(bounds.inside({0, 2}, tolerance));
    BOOST_CHECK(bounds.inside({2, 0}, tolerance));
  }

  {
    auto tolerance = BoundaryTolerance::Chi2Bound(SquareMatrix2::Identity(), 2);
    BOOST_CHECK(bounds.inside({0, 0}, tolerance));
    BOOST_CHECK(bounds.inside({2, 2}, tolerance));
    BOOST_CHECK(bounds.inside({0, 2}, tolerance));
    BOOST_CHECK(bounds.inside({2, 0}, tolerance));
    BOOST_CHECK(!bounds.inside({3, 3}, tolerance));
  }
}

BOOST_AUTO_TEST_CASE(BoundaryCheckNegativeToleranceRect) {
  // Test points for boundary check with euclidean tolerance
  Vector2 ll(1, 1);
  Vector2 ur(3, 3);
  RectangleBounds bounds(ll, ur);

  auto check = [&bounds](const BoundaryTolerance& tolerance,
                         const Vector2& point) {
    return bounds.inside(point, tolerance);
  };

  {
    auto tolerance = BoundaryTolerance::AbsoluteEuclidean(-0.25);

    BOOST_CHECK(!check(tolerance, {2.8, 2}));
    BOOST_CHECK(!check(tolerance, {3.1, 2}));
    BOOST_CHECK(check(tolerance, {2.7, 2}));
    BOOST_CHECK(!check(tolerance, {2, 3.1}));
    BOOST_CHECK(!check(tolerance, {2, 2.8}));
    BOOST_CHECK(check(tolerance, {2, 2.7}));

    BOOST_CHECK(!check(tolerance, {0.8, 2}));
    BOOST_CHECK(!check(tolerance, {1.2, 2}));
    BOOST_CHECK(check(tolerance, {1.5, 2}));
    BOOST_CHECK(!check(tolerance, {2, 0.8}));
    BOOST_CHECK(!check(tolerance, {2, 1.2}));
    BOOST_CHECK(check(tolerance, {2, 1.5}));
  }

  {
    auto tolerance =
        BoundaryTolerance::Chi2Bound(SquareMatrix2::Identity(), -0.2);

    BOOST_CHECK(!check(tolerance, {2.8, 2}));
    BOOST_CHECK(!check(tolerance, {3.1, 2}));
    BOOST_CHECK(check(tolerance, {2.5, 2}));
    BOOST_CHECK(!check(tolerance, {2, 3.1}));
    BOOST_CHECK(!check(tolerance, {2, 2.8}));
    BOOST_CHECK(check(tolerance, {2, 2.5}));

    BOOST_CHECK(!check(tolerance, {0.8, 2}));
    BOOST_CHECK(!check(tolerance, {1.4, 2}));
    BOOST_CHECK(check(tolerance, {1.5, 2}));
    BOOST_CHECK(!check(tolerance, {2, 0.8}));
    BOOST_CHECK(!check(tolerance, {2, 1.4}));
    BOOST_CHECK(check(tolerance, {2, 1.5}));
  }
}

BOOST_AUTO_TEST_CASE(BoundaryCheckNegativeToleranceTrap) {
  ConvexPolygonBounds<PolygonDynamic> bounds(
      {{1.5, 1}, {2.5, 1}, {3, 3}, {1, 3}});

  auto check = [&bounds](const BoundaryTolerance& tolerance,
                         const Vector2& point) {
    return bounds.inside(point, tolerance);
  };

  {
    auto tolerance = BoundaryTolerance::AbsoluteEuclidean(0.25);
    // Axes
    BOOST_CHECK(!check(tolerance, {3.1, 2}));
    BOOST_CHECK(check(tolerance, {2.75, 2}));
    BOOST_CHECK(check(tolerance, {2.5, 2}));
    BOOST_CHECK(check(tolerance, {2.25, 2}));
    BOOST_CHECK(check(tolerance, {2, 3.1}));
    BOOST_CHECK(check(tolerance, {2, 2.75}));
    BOOST_CHECK(check(tolerance, {2, 2.5}));
    BOOST_CHECK(check(tolerance, {2, 2.25}));
    BOOST_CHECK(check(tolerance, {2, 2}));

    // Corners
    BOOST_CHECK(check(tolerance, {3.1, 3.2}));
    BOOST_CHECK(check(tolerance, {0.9, 3.2}));
    BOOST_CHECK(check(tolerance, {1.5, 0.8}));
    BOOST_CHECK(check(tolerance, {2.5, 0.8}));
  }

  {
    auto tolerance = BoundaryTolerance::AbsoluteEuclidean(-0.25);
    // Axes
    BOOST_CHECK(!check(tolerance, {3.0, 2}));
    BOOST_CHECK(!check(tolerance, {2.5, 2}));
    BOOST_CHECK(check(tolerance, {2.25, 2}));
    BOOST_CHECK(!check(tolerance, {2, 3.1}));
    BOOST_CHECK(!check(tolerance, {2, 2.9}));
    BOOST_CHECK(check(tolerance, {2, 2.7}));

    // Corners
    BOOST_CHECK(!check(tolerance, {2.7, 2.9}));
    BOOST_CHECK(check(tolerance, {2.4, 2.6}));
    BOOST_CHECK(!check(tolerance, {1.3, 2.9}));
    BOOST_CHECK(check(tolerance, {1.6, 2.6}));
    BOOST_CHECK(!check(tolerance, {2.4, 1.1}));
    BOOST_CHECK(check(tolerance, {1.75, 1.4}));
    BOOST_CHECK(!check(tolerance, {1.6, 1.1}));
    BOOST_CHECK(check(tolerance, {2.25, 1.4}));
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
