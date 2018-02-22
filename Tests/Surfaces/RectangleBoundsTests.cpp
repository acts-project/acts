// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Rectangle Bounds Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include <algorithm>
#include <limits>
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace utf    = boost::unit_test;
const double inf = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  void
  dumpVertices(const RectangleBounds& r)
  {
    const auto& v = r.vertices();
    for (const auto& i : v) {
      std::cout << "(" << i[0] << ", " << i[1] << ")" << std::endl;
    }
  }
  bool
  approximatelyEqual(const Vector2D& a, const Vector2D& b)
  {
    const double dif0 = std::abs(a[0] - b[0]);
    const double dif1 = std::abs(a[1] - b[1]);
    const double tol  = 1e-9;
    return ((dif0 < tol) and (dif1 < tol));
  }
  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant RectangleBounds object
  BOOST_AUTO_TEST_CASE(RectangleBoundsConstruction)
  {
    const double    halfX(10.), halfY(5.);
    RectangleBounds twentyByTenRectangle(halfX, halfY);
    BOOST_TEST(twentyByTenRectangle.type() == Acts::SurfaceBounds::Rectangle);
    //
    // nonsensical bounds are also permitted, but maybe should not be
    const double zeroHalfX(0.), zeroHalfY(0.);
    const double infHalfX(inf), infHalfY(inf);
    const double nanHalfX(NaN), nanHalfY(NaN);
    const double negHalfX(-10.), negHalfY(-5.);
    //
    // BOOST_TEST_MESSAGE("Initialise with zero dimensions");
    RectangleBounds zeroDimensionsRectangle(zeroHalfX, zeroHalfY);
    BOOST_TEST(zeroDimensionsRectangle.type()
               == Acts::SurfaceBounds::Rectangle);
    //
    // BOOST_TEST_MESSAGE("Initialise with infinite dimensions");
    RectangleBounds infinite(infHalfX, infHalfY);
    BOOST_TEST(infinite.type() == Acts::SurfaceBounds::Rectangle);
    //
    // BOOST_TEST_MESSAGE("Initialise with NaN dimensions");
    RectangleBounds nanRectangle(nanHalfX, nanHalfY);
    BOOST_TEST(nanRectangle.type() == Acts::SurfaceBounds::Rectangle);
    //
    // BOOST_TEST_MESSAGE("Initialise with negative dimensions");
    RectangleBounds negativeDimensionedRectangle(negHalfX, negHalfY);
    BOOST_TEST(negativeDimensionedRectangle.type()
               == Acts::SurfaceBounds::Rectangle);
  }

  /// Unit test for testing RectangleBounds properties
  BOOST_TEST_DECORATOR(*utf::tolerance(1e-10))
  BOOST_AUTO_TEST_CASE(RectangleBoundsProperties)
  {
    const double    halfX(10.), halfY(5.);
    RectangleBounds rect(halfX, halfY);
    BOOST_TEST(rect.halflengthX() == 10.);
    BOOST_TEST(rect.halflengthY() == 5.);
    const std::vector<Vector2D> coords
        = {{10., -5.}, {10., 5.}, {-10., 5.}, {-10., -5.}};
    // equality, ensure ordering is ok
    BOOST_TEST(
        std::equal(coords.begin(), coords.end(), rect.vertices().begin()));
    const Vector2D pointA{1.0, 1.0}, pointB{9.0, 1.0}, outside{10.1, 5.1};
    // distance is signed, from boundary to point. (doesn't seem right, given
    // the name of the method)
    BOOST_TEST(rect.distanceToBoundary(pointA) == -4.0);
    BOOST_TEST(rect.distanceToBoundary(pointB) == -1.0);
    BoundaryCheck bcheck(true, true);
    BOOST_TEST(rect.inside(pointA, bcheck));
  }
  BOOST_AUTO_TEST_CASE(RectangleBoundsAssignment)
  {
    const double    halfX(10.), halfY(2.);
    RectangleBounds rectA(halfX, halfY);
    RectangleBounds rectB(0.0, 0.0);
    rectB                       = rectA;
    const auto originalVertices = rectA.vertices();
    const auto assignedVertices = rectB.vertices();
    BOOST_TEST(originalVertices == assignedVertices);
  }

  BOOST_AUTO_TEST_CASE(RectangleBoundsClone)
  {
    const double    halfX(10.), halfY(5.);
    RectangleBounds rectA(halfX, halfY);
    auto            rectB = rectA.clone();
    BOOST_TEST(bool(rectB));  // not null pointer
    const auto& originalVertices = rectA.vertices();
    const auto& clonedVertices   = rectB->vertices();
    BOOST_TEST(originalVertices == clonedVertices);
    delete rectB;
  }
  BOOST_AUTO_TEST_SUITE_END()
}  // end of namespace Test

}  // end of namespace Acts
