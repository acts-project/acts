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

#include <limits>

#include "Acts/Surfaces/PlanarBooleanBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit tests for PlanarBooleanBounds - union type
BOOST_AUTO_TEST_CASE(PlanarBooleanBoundsUnion) {
  auto lBounds = std::make_shared<const RectangleBounds>(3, 2.5);
  auto rBounds = std::make_shared<const RectangleBounds>(1.5, 1.);

  PlanarBooleanBounds recSepUnion(lBounds, Vector2D(0., -4), rBounds,
                                  Vector2D(0., 3.), eUnion);

  // Check point inside the first
  BOOST_CHECK(recSepUnion.inside({-2., -3.}, true));
  // Check point inside the second
  BOOST_CHECK(recSepUnion.inside({0.5, 3.5}, true));
  // Check point outside both
  BOOST_CHECK(!recSepUnion.inside({0.5, 0.}, true));
  // Check the overlap
  BOOST_CHECK(!recSepUnion.boundsOverlap());
  // Check the bounding box
  auto recSepUnionBB = recSepUnion.boundingBox();
  auto& recSepUnionBBmin = recSepUnionBB.min();
  auto& recSepUnionBBmax = recSepUnionBB.max();
  CHECK_CLOSE_ABS(recSepUnionBBmin.x(), -3., s_epsilon);
  CHECK_CLOSE_ABS(recSepUnionBBmin.y(), -6.5, s_epsilon);
  CHECK_CLOSE_ABS(recSepUnionBBmax.x(), 3., s_epsilon);
  CHECK_CLOSE_ABS(recSepUnionBBmax.y(), 4., s_epsilon);

  // Check the distance to boundary - 1.5 away from left
  CHECK_CLOSE_ABS(recSepUnion.distanceToBoundary({0.5, 0.}), 1.5, s_epsilon);

  PlanarBooleanBounds recOvUnion(lBounds, Vector2D(0., 0.), rBounds,
                                 Vector2D(3., 0.), eUnion);

  // Check the overlap
  BOOST_CHECK(recOvUnion.boundsOverlap());
  // Check point inside the intersection
  BOOST_CHECK(recOvUnion.inside({1.5, 0.5}, true));
  // Check point outside the intersection, but inside first
  BOOST_CHECK(recOvUnion.inside({0., 0.}, true));
  // Check point outside the intersection, but inside second
  BOOST_CHECK(recOvUnion.inside({4., 0.}, true));
}

/// Unit tests for PlanarBooleanBounds - intersection type
BOOST_AUTO_TEST_CASE(PlanarBooleanBoundsIntersection) {
  auto lBounds = std::make_shared<const RectangleBounds>(3, 2.5);
  auto rBounds = std::make_shared<const RectangleBounds>(2., 1.);

  PlanarBooleanBounds recIntersection(lBounds, Vector2D(0., 0.), rBounds,
                                      Vector2D(3., 0.), eIntersection);

  // Check the overlap
  BOOST_CHECK(recIntersection.boundsOverlap());
  // Check point inside the intersection
  BOOST_CHECK(recIntersection.inside({1.5, 0.5}, true));
  // Check point outside the intersection, but inside first
  BOOST_CHECK(!recIntersection.inside({0., 0.}, true));
  // Check point outside the intersection, but inside second
  BOOST_CHECK(!recIntersection.inside({5., 0.}, true));

  // Check the distance to boundary - 1.5 away from left
  // CHECK_CLOSE_ABS(recIntersection.distanceToBoundary({0., 0.}), -1.,
  // s_epsilon);

  BOOST_CHECK_THROW(PlanarBooleanBounds(lBounds, Vector2D(0., 0.), rBounds,
                                        Vector2D(10., 0.), eIntersection),
                    std::logic_error);
}

/// Unit tests for PlanarBooleanBounds - not type
BOOST_AUTO_TEST_CASE(PlanarBooleanBoundsNot) {
  auto lBounds = std::make_shared<const RectangleBounds>(3, 2.5);
  auto rBounds = std::make_shared<const RectangleBounds>(2., 1.);

  PlanarBooleanBounds recNot(lBounds, Vector2D(0., 0.), rBounds,
                             Vector2D(3., 0.), eNot);

  // Check the overlap
  BOOST_CHECK(recNot.boundsOverlap());
  // Check point inside the cutout
  BOOST_CHECK(!recNot.inside({1.5, 0.5}, true));
  // Check point outside the cutout, but inside first
  BOOST_CHECK(recNot.inside({0., 0.}, true));
  // Check point outside the cut, but inside second
  BOOST_CHECK(!recNot.inside({5., 0.}, true));

  // Check the distance to boundary - 1.5 away from left
  // CHECK_CLOSE_ABS(recIntersection.distanceToBoundary({0., 0.}), -1.,
  // s_epsilon);

  BOOST_CHECK_THROW(PlanarBooleanBounds(lBounds, Vector2D(0., 0.), rBounds,
                                        Vector2D(10., 0.), eNot),
                    std::logic_error);
}

/// Unit tests for PlanarBooleanBounds - XOR type
BOOST_AUTO_TEST_CASE(PlanarBooleanBoundsXor) {
  auto lBounds = std::make_shared<const RectangleBounds>(3, 2.5);
  auto rBounds = std::make_shared<const RectangleBounds>(2., 1.);

  PlanarBooleanBounds recXor(lBounds, Vector2D(0., 0.), rBounds,
                             Vector2D(3., 0.), eXor);

  // Check the overlap
  BOOST_CHECK(recXor.boundsOverlap());
  // Check point inside the cutout
  BOOST_CHECK(!recXor.inside({1.5, 0.5}, true));
  // Check point outside the cutout, but inside first
  BOOST_CHECK(recXor.inside({0., 0.}, true));
  // Check point outside the cut, but inside second
  BOOST_CHECK(recXor.inside({4., 0.}, true));

  // Check the distance to boundary - 1.5 away from left
  // CHECK_CLOSE_ABS(recIntersection.distanceToBoundary({0., 0.}), -1.,
  // s_epsilon);

  BOOST_CHECK_THROW(PlanarBooleanBounds(lBounds, Vector2D(0., 0.), rBounds,
                                        Vector2D(10., 0.), eXor),
                    std::logic_error);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
