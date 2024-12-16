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
#include "Acts/Surfaces/detail/BoundaryCheckHelper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <limits>
#include <optional>
#include <random>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Surfaces)

// See: https://en.wikipedia.org/wiki/Bounding_volume
//
// Aligned box w/ simple check
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxSimple) {
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  auto tolerance = BoundaryTolerance::None();
  BOOST_CHECK(
      detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
  BOOST_CHECK(
      !detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
  BOOST_CHECK(
      !detail::insideAlignedBox(ll, ur, tolerance, {0, 2}, std::nullopt));
  BOOST_CHECK(
      !detail::insideAlignedBox(ll, ur, tolerance, {2, 0}, std::nullopt));
}

// Aligned box w/ tolerance check along first axis
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxToleranceLoc0) {
  boost::execution_monitor em;
  em.p_detect_fp_exceptions.set(boost::fpe::BOOST_FPE_ALL);

  em.execute([]() {
    Vector2 ll(-1, -1);
    Vector2 ur(1, 1);
    auto tolerance = BoundaryTolerance::AbsoluteBound(
        1.5, std::numeric_limits<double>::infinity());
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {4, 4}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 2}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {2, 0}, std::nullopt));

    return 0;
  });
}

// Aligned box w/ covariance check
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxCovariance) {
  SquareMatrix2 cov;
  cov << 1, 0.5, 0.5, 2;
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  auto tolerance = BoundaryTolerance::Chi2Bound(cov.inverse(), 3.);
  BOOST_CHECK(
      detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
  BOOST_CHECK(
      detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
  BOOST_CHECK(
      !detail::insideAlignedBox(ll, ur, tolerance, {4, 4}, std::nullopt));
  BOOST_CHECK(
      detail::insideAlignedBox(ll, ur, tolerance, {0, 3}, std::nullopt));
  BOOST_CHECK(
      detail::insideAlignedBox(ll, ur, tolerance, {3, 0}, std::nullopt));
}

// Triangle w/ simple check
BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleSimple) {
  Vector2 vertices[] = {{-2, 0}, {2, 0}, {0, 2}};
  auto tolerance = BoundaryTolerance::None();
  BOOST_CHECK(detail::insidePolygon(vertices, tolerance, {0, 0}, std::nullopt));
  BOOST_CHECK(detail::insidePolygon(vertices, tolerance, {0, 1}, std::nullopt));
  BOOST_CHECK(
      !detail::insidePolygon(vertices, tolerance, {2, 2}, std::nullopt));
  BOOST_CHECK(
      !detail::insidePolygon(vertices, tolerance, {0, -1}, std::nullopt));
}
// Triangle w/ covariance check
BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleCovariance) {
  Vector2 vertices[] = {{-2, 0}, {2, 0}, {0, 2}};
  SquareMatrix2 cov;
  cov << 0.5, 0, 0, 0.5;
  auto tolerance = BoundaryTolerance::Chi2Bound(cov.inverse(), 4.1);
  BOOST_CHECK(detail::insidePolygon(vertices, tolerance, {0, 0}, std::nullopt));
  BOOST_CHECK(detail::insidePolygon(vertices, tolerance, {0, 1}, std::nullopt));
  BOOST_CHECK(detail::insidePolygon(vertices, tolerance, {0, 2}, std::nullopt));
  BOOST_CHECK(detail::insidePolygon(vertices, tolerance, {0, 3}, std::nullopt));
  BOOST_CHECK(detail::insidePolygon(vertices, tolerance, {0, 4}, std::nullopt));
  BOOST_CHECK(
      !detail::insidePolygon(vertices, tolerance, {0, 5}, std::nullopt));
}

BOOST_AUTO_TEST_CASE(BoundaryCheckDifferentTolerances) {
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);

  {
    auto tolerance = BoundaryTolerance::None();
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {0, 2}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {2, 0}, std::nullopt));
  }

  {
    auto tolerance = BoundaryTolerance::Infinite();
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 2}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {2, 0}, std::nullopt));
  }

  {
    auto tolerance = BoundaryTolerance::AbsoluteBound(0.5, 0.5);
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {1.1, 1.1}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 1.1}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {1.1, 0}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {0, 2}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {2, 0}, std::nullopt));
  }

  {
    auto tolerance = BoundaryTolerance::AbsoluteCartesian(0.5, 0.5);
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {1.1, 1.1}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 1.1}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {1.1, 0}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {0, 2}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {2, 0}, std::nullopt));
  }

  {
    auto tolerance = BoundaryTolerance::AbsoluteEuclidean(1.1);
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 2}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {2, 0}, std::nullopt));
  }

  {
    auto tolerance =
        BoundaryTolerance::Chi2Bound(SquareMatrix2::Identity(), 1.);
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {0, 2}, std::nullopt));
    BOOST_CHECK(
        detail::insideAlignedBox(ll, ur, tolerance, {2, 0}, std::nullopt));
    BOOST_CHECK(
        !detail::insideAlignedBox(ll, ur, tolerance, {3, 3}, std::nullopt));
  }
}

BOOST_AUTO_TEST_CASE(BoundaryCheckNegativeTolerance) {
  // Test points for boundary check with euclidean tolerance
  Vector2 ll(1, 1);
  Vector2 ur(3, 3);

  auto check = [&ll, &ur](const BoundaryTolerance& tolerance,
                          const Vector2& point) {
    return detail::insideAlignedBox(ll, ur, tolerance, point, std::nullopt);
  };

  {
    auto tolerance = BoundaryTolerance::AbsoluteEuclidean(-0.25);

    BOOST_CHECK(!check(tolerance, {2.8, 2}));
    BOOST_CHECK(!check(tolerance, {3.1, 2}));
    BOOST_CHECK(check(tolerance, {2.7, 2}));
    BOOST_CHECK(!check(tolerance, {2, 3.1}));
    BOOST_CHECK(!check(tolerance, {2, 2.8}));
    BOOST_CHECK(check(tolerance, {2, 2.7}));
  }

  {
    auto tolerance = BoundaryTolerance::AbsoluteBound(-0.25, -0.0);

    BOOST_CHECK(!check(tolerance, {2.8, 2}));
    BOOST_CHECK(!check(tolerance, {3.1, 2}));
    BOOST_CHECK(check(tolerance, {2.7, 2}));
    BOOST_CHECK(!check(tolerance, {2, 3.1}));
    BOOST_CHECK(check(tolerance, {2, 2.8}));
    BOOST_CHECK(check(tolerance, {2, 2.7}));
  }

  // std::ofstream outFile("boundary_check_euclidean.csv");
  // outFile << "x,y,inone,ipos,ineg\n";
  //
  // std::mt19937 gen(42);
  // std::uniform_real_distribution<> dis(0.5, 3.5);
  //
  // for (int i = 0; i < 1000; i++) {
  //   double x = dis(gen);
  //   double y = dis(gen);
  //   Vector2 point(x, y);
  //
  //   bool insideNone =
  //       detail::insideAlignedBox(ll, ur, none, point, std::nullopt);
  //   bool insidePos = detail::insideAlignedBox(ll, ur, pos, point,
  //   std::nullopt); bool insideNeg = detail::insideAlignedBox(ll, ur, neg,
  //   point, std::nullopt);
  //
  //   outFile << x << "," << y << "," << insideNone << "," << insidePos << ","
  //           << insideNeg << "\n";
  // }
  //
  // outFile.close();
  //
  // return;
  // Verify some key points with known distances
  // BOOST_CHECK(
  //     detail::insideAlignedBox(ll, ur, tolerance, {0, 0}, std::nullopt));
  // BOOST_CHECK(detail::insideAlignedBox(ll, ur, tolerance, {1.05, 0},
  //                                      std::nullopt));  // Just within
  //                                      tolerance
  // BOOST_CHECK(!detail::insideAlignedBox(
  //     ll, ur, tolerance, {1.2, 0}, std::nullopt));  // Just outside tolerance
  // BOOST_CHECK(
  //     !detail::insideAlignedBox(ll, ur, tolerance, {2, 2}, std::nullopt));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
