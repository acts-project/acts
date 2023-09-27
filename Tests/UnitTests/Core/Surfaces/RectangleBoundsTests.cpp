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

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <iostream>
#include <limits>

namespace utf = boost::unit_test;
const double inf = std::numeric_limits<double>::infinity();

namespace Acts {

namespace Test {
void dumpVertices(const RectangleBounds& r) {
  const auto& v = r.vertices();
  for (const auto& i : v) {
    std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
  }
}
bool approximatelyEqual(const Vector2& a, const Vector2& b) {
  const double dif0 = std::abs(a[0] - b[0]);
  const double dif1 = std::abs(a[1] - b[1]);
  const double tol = 1e-9;
  return ((dif0 < tol) and (dif1 < tol));
}
BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit test for creating compliant/non-compliant RectangleBounds object
BOOST_AUTO_TEST_CASE(RectangleBoundsConstruction) {
  const double halfX(10.), halfY(5.);
  RectangleBounds twentyByTenRectangle(halfX, halfY);
  BOOST_CHECK_EQUAL(twentyByTenRectangle.type(),
                    Acts::SurfaceBounds::eRectangle);
  //
  // nonsensical bounds are also permitted, but maybe should not be
  const double zeroHalfX(0.), zeroHalfY(0.);
  const double infHalfX(inf), infHalfY(inf);
  //
  // BOOST_TEST_MESSAGE("Initialise with zero dimensions");
  RectangleBounds zeroDimensionsRectangle(zeroHalfX, zeroHalfY);
  BOOST_CHECK_EQUAL(zeroDimensionsRectangle.type(),
                    Acts::SurfaceBounds::eRectangle);
  //
  // BOOST_TEST_MESSAGE("Initialise with infinite dimensions");
  RectangleBounds infinite(infHalfX, infHalfY);
  BOOST_CHECK_EQUAL(infinite.type(), Acts::SurfaceBounds::eRectangle);
}

/// Recreation
BOOST_AUTO_TEST_CASE(RectangleBoundsRecreation) {
  const double halfX(10.), halfY(2.);
  RectangleBounds original(halfX, halfY);
  // const bool symmetric(false);
  auto valvector = original.values();
  std::array<double, RectangleBounds::eSize> values{};
  std::copy_n(valvector.begin(), RectangleBounds::eSize, values.begin());
  RectangleBounds recreated(values);
  BOOST_CHECK_EQUAL(original, recreated);
}

// Exception tests
BOOST_AUTO_TEST_CASE(RadialBoundsException) {
  const double halfX(10.), halfY(2.);

  // Negative x half length
  BOOST_CHECK_THROW(RectangleBounds(-halfX, halfY), std::logic_error);

  // Negative y half length
  BOOST_CHECK_THROW(RectangleBounds(halfX, -halfY), std::logic_error);
}

/// Unit test for testing RectangleBounds properties
BOOST_TEST_DECORATOR(*utf::tolerance(1e-10))
BOOST_AUTO_TEST_CASE(RectangleBoundsProperties) {
  const double halfX(10.), halfY(5.);
  RectangleBounds rect(halfX, halfY);
  BOOST_CHECK_EQUAL(rect.halfLengthX(), 10.);
  BOOST_CHECK_EQUAL(rect.halfLengthY(), 5.);

  CHECK_CLOSE_ABS(rect.min(), Vector2(-halfX, -halfY), 1e-6);
  CHECK_CLOSE_ABS(rect.max(), Vector2(halfX, halfY), 1e-6);

  const std::vector<Vector2> coords = {
      {-10., -5.}, {10., -5.}, {10., 5.}, {-10., 5.}};
  // equality, ensure ordering is ok
  const auto& rectVertices = rect.vertices();
  BOOST_CHECK_EQUAL_COLLECTIONS(coords.cbegin(), coords.cend(),
                                rectVertices.cbegin(), rectVertices.cend());
  const Vector2 pointA{1.0, 1.0};
  // distance is signed, from boundary to point. (doesn't seem right, given
  BoundaryCheck bcheck(true, true);
  BOOST_CHECK(rect.inside(pointA, bcheck));
}
BOOST_AUTO_TEST_CASE(RectangleBoundsAssignment) {
  const double halfX(10.), halfY(2.);
  RectangleBounds rectA(halfX, halfY);
  RectangleBounds rectB(0.0, 0.0);
  rectB = rectA;
  const auto originalVertices = rectA.vertices();
  const auto assignedVertices = rectB.vertices();
  BOOST_CHECK_EQUAL_COLLECTIONS(
      originalVertices.cbegin(), originalVertices.cend(),
      assignedVertices.cbegin(), assignedVertices.cend());
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test

}  // namespace Acts
