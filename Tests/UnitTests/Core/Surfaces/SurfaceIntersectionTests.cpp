// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <cmath>
#include <memory>
#include <utility>

namespace Acts {

using namespace UnitLiterals;

// Create a test context
GeometryContext tgContext = GeometryContext();

// Some random transform
Transform3 aTransform = Transform3::Identity() *
                        Translation3(30_cm, 7_m, -87_mm) *
                        AngleAxis3(0.42, Vector3(-3., 1., 8).normalized());

namespace Test {

BOOST_AUTO_TEST_SUITE(Surfaces)

/// This tests the intersection with cylinders
/// and looks for valid, non-valid, solutions
BOOST_AUTO_TEST_CASE(CylinderIntersectionTests) {
  double radius = 1_m;
  double halfZ = 10_m;

  auto testCylinderIntersection = [&](const Transform3& transform) -> void {
    // A cylinder created aligned with a provided transform
    auto aCylinder =
        Surface::makeShared<CylinderSurface>(transform, radius, halfZ);

    // Linear transform
    auto lTransform = transform.linear();

    // An onCylinder solution
    Vector3 onCylinder = transform * Vector3(radius, 0., 0.);
    Vector3 outCylinder = transform * Vector3(-radius, 0.6 * radius, 90_cm);
    Vector3 atCenter = transform * Vector3(0., 0., 0.);
    Vector3 atEdge = transform * Vector3(0.5 * radius, 0., 0.99 * halfZ);
    // Simply along the x axis
    Vector3 alongX = lTransform * Vector3(1., 0., 0.);
    Vector3 transXY = lTransform * Vector3(1., 1., 0).normalized();
    Vector3 transTZ = lTransform * Vector3(1., 0., 1.).normalized();

    // Intersect without boundary check
    auto aIntersection =
        aCylinder->intersect(tgContext, onCylinder, alongX, true);

    // Check the validity of the intersection
    BOOST_CHECK(aIntersection[0]);
    // The status of this one should be on surface
    BOOST_CHECK(aIntersection[0].status() == Intersection3D::Status::reachable);
    // The intersection is at 2 meter distance
    CHECK_CLOSE_ABS(aIntersection[0].pathLength(), -2_m, s_onSurfaceTolerance);
    // There MUST be a second solution
    BOOST_CHECK(aIntersection[1]);
    // The other intersection MUST be reachable
    BOOST_CHECK(aIntersection[1].status() == Intersection3D::Status::onSurface);

    // Intersect from the center
    auto cIntersection =
        aCylinder->intersect(tgContext, atCenter, alongX, true);

    // Check the validity of the intersection
    BOOST_CHECK(cIntersection[0]);
    // The status of this one MUST be reachable
    BOOST_CHECK(cIntersection[0].status() == Intersection3D::Status::reachable);
    // There MUST be a second solution
    BOOST_CHECK(cIntersection[1]);
    // The other intersection MUST be reachable
    BOOST_CHECK(cIntersection[1].status() == Intersection3D::Status::reachable);
    // There MUST be one forward one backwards solution
    BOOST_CHECK(cIntersection[1].pathLength() * cIntersection[0].pathLength() <
                0);

    // Intersect from outside where both intersections are reachable
    auto oIntersection =
        aCylinder->intersect(tgContext, outCylinder, alongX, true);

    // Check the validity of the intersection
    BOOST_CHECK(oIntersection[0]);
    // The status of this one MUST be reachable
    BOOST_CHECK(oIntersection[0].status() == Intersection3D::Status::reachable);
    // There MUST be a second solution
    BOOST_CHECK(oIntersection[1]);
    // The other intersection MUST be reachable
    BOOST_CHECK(oIntersection[1].status() == Intersection3D::Status::reachable);
    // There MUST be one forward one backwards solution
    BOOST_CHECK(oIntersection[1].pathLength() * oIntersection[0].pathLength() >
                0);

    // Intersection from outside without chance of hitting the cylinder
    auto iIntersection =
        aCylinder->intersect(tgContext, outCylinder, transXY, false);

    // Check the validity of the intersection
    BOOST_CHECK(!iIntersection[0]);

    // From edge tests - wo boundary test
    auto eIntersection =
        aCylinder->intersect(tgContext, atEdge, transTZ, false);

    // Check the validity of the intersection
    BOOST_CHECK(eIntersection[0]);
    // This should be the positive one
    BOOST_CHECK(eIntersection[0].pathLength() < 0.);
    // The status of this one should be reachable
    BOOST_CHECK(eIntersection[0].status() == Intersection3D::Status::reachable);
    // There MUST be a second solution
    BOOST_CHECK(eIntersection[1]);
    // The other intersection MUST be reachable
    BOOST_CHECK(eIntersection[1].status() == Intersection3D::Status::reachable);
    // And be the negative one
    BOOST_CHECK(eIntersection[1].pathLength() > 0.);

    // Now re-do with boundary check
    eIntersection = aCylinder->intersect(tgContext, atEdge, transTZ, true);
    // This should be the negative one
    BOOST_CHECK(eIntersection[0].pathLength() < 0.);
    // The status of this one should be reachable
    BOOST_CHECK(eIntersection[0].status() == Intersection3D::Status::reachable);
    // There MUST be a second solution
    BOOST_CHECK(!eIntersection[1]);
    // The other intersection MUST NOT be reachable
    BOOST_CHECK(eIntersection[1].status() == Intersection3D::Status::missed);
    // And be the positive one
    BOOST_CHECK(eIntersection[1].pathLength() > 0.);
  };

  // In a nominal world
  testCylinderIntersection(Transform3::Identity());

  // In a system somewhere away
  testCylinderIntersection(aTransform);
}

/// This tests the intersection with cylinders
/// and looks for valid, non-valid, solutions
BOOST_AUTO_TEST_CASE(ConeIntersectionTest) {
  double alpha = 0.25 * M_PI;

  auto testConeIntersection = [&](const Transform3& transform) -> void {
    // A cone surface ready to use
    auto aCone = Surface::makeShared<ConeSurface>(transform, alpha, true);

    // Linear transform
    auto lTransform = transform.linear();

    // An onCylinder solution
    Vector3 onCone = transform * Vector3(std::sqrt(2.), std::sqrt(2.), 2.);
    Vector3 outCone = transform * Vector3(std::sqrt(4.), std::sqrt(4.), 2.);
    // Simply along the x axis
    Vector3 perpXY = lTransform * Vector3(1., -1., 0.).normalized();
    Vector3 transXY = lTransform * Vector3(1., 1., 0).normalized();

    // Intersect without boundary check with an on solution
    BOOST_CHECK(aCone->isOnSurface(tgContext, onCone, transXY, false));
    auto aIntersection = aCone->intersect(tgContext, onCone, transXY, true);

    // Check the validity of the intersection
    BOOST_CHECK(aIntersection[0]);
    // The status of this one should be on surface
    BOOST_CHECK(aIntersection[0].status() == Intersection3D::Status::reachable);
    // The intersection is at 4 mm distance
    CHECK_CLOSE_ABS(aIntersection[0].pathLength(), -4., s_onSurfaceTolerance);
    // There MUST be a second solution
    BOOST_CHECK(aIntersection[1]);
    // The other intersection MUST be reachable
    BOOST_CHECK(aIntersection[1].status() == Intersection3D::Status::onSurface);

    // Intersection from outside without chance of hitting the cylinder
    auto iIntersection = aCone->intersect(tgContext, outCone, perpXY, false);

    // Check the validity of the intersection
    BOOST_CHECK(!iIntersection[0]);
  };

  // In a nominal world
  testConeIntersection(Transform3::Identity());

  // In a system somewhere away
  testConeIntersection(aTransform);
}

/// This tests the intersection with planar surfaces (plane, disk)
/// as those share the same PlanarHelper methods, only one test is
/// sufficient
/// - it looks for valid, non-valid, solutions
BOOST_AUTO_TEST_CASE(PlanarIntersectionTest) {
  double halfX = 1_m;
  double halfY = 10_m;

  auto testPlanarIntersection = [&](const Transform3& transform) -> void {
    // A Plane created with a specific transform
    auto aPlane = Surface::makeShared<PlaneSurface>(
        transform, std::make_shared<RectangleBounds>(halfX, halfY));

    /// Forward intersection test
    Vector3 before = transform * Vector3(-50_cm, -1_m, -1_m);
    Vector3 onit = transform * Vector3(11_cm, -22_cm, 0_m);
    Vector3 after = transform * Vector3(33_cm, 12_mm, 1_m);
    Vector3 outside = transform * Vector3(2. * halfX, 2 * halfY, -1_mm);

    // Linear transform
    auto lTransform = transform.linear();

    // A direction that is non trivial
    Vector3 direction = lTransform * Vector3(4_mm, 8_mm, 50_cm).normalized();
    Vector3 parallel = lTransform * Vector3(1., 1., 0.).normalized();

    // Intersect forward
    auto fIntersection = aPlane->intersect(tgContext, before, direction, true);

    // The intersection MUST be valid
    BOOST_CHECK(fIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(fIntersection[0].status() == Intersection3D::Status::reachable);
    // The path length MUST be positive
    BOOST_CHECK(fIntersection[0].pathLength() > 0.);
    // The intersection MUST be unique
    BOOST_CHECK(!fIntersection[1]);

    // On surface intersection
    auto oIntersection = aPlane->intersect(tgContext, onit, direction, true);
    // The intersection MUST be valid
    BOOST_CHECK(oIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(oIntersection[0].status() == Intersection3D::Status::onSurface);
    // The path length MUST be positive
    BOOST_CHECK(std::abs(oIntersection[0].pathLength()) < s_onSurfaceTolerance);
    // The intersection MUST be unique
    BOOST_CHECK(!oIntersection[1]);

    // Intersect backwards
    auto bIntersection = aPlane->intersect(tgContext, after, direction, true);
    // The intersection MUST be valid
    BOOST_CHECK(bIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(bIntersection[0].status() == Intersection3D::Status::reachable);
    // The path length MUST be negative
    BOOST_CHECK(bIntersection[0].pathLength() < 0.);
    // The intersection MUST be unique
    BOOST_CHECK(!bIntersection[1]);

    // An out of bounds attempt: missed
    auto mIntersection = aPlane->intersect(tgContext, outside, direction, true);
    // The intersection MUST NOT be valid
    BOOST_CHECK(!mIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(mIntersection[0].status() == Intersection3D::Status::missed);
    // The path length MUST be negative
    BOOST_CHECK(mIntersection[0].pathLength() > 0.);
    // The intersection MUST be unique
    BOOST_CHECK(!mIntersection[1]);

    // An invalid attempt
    auto iIntersection = aPlane->intersect(tgContext, before, parallel, true);
    // The intersection MUST NOT be valid
    BOOST_CHECK(!iIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(iIntersection[0].status() ==
                Intersection3D::Status::unreachable);
    // The intersection MUST be unique
    BOOST_CHECK(!iIntersection[1]);
  };

  // In a nominal world
  testPlanarIntersection(Transform3::Identity());

  // In a system somewhere away
  testPlanarIntersection(aTransform);
}

/// This tests the intersection with line like surfaces (straw, perigee)
/// as those share the same methods, only one test is
/// sufficient
/// - it looks for valid, non-valid, solutions
BOOST_AUTO_TEST_CASE(LineIntersectionTest) {
  double radius = 1_m;
  double halfZ = 10_m;

  auto testLineAppraoch = [&](const Transform3& transform) -> void {
    // A Plane created with a specific transform
    auto aLine = Surface::makeShared<StrawSurface>(transform, radius, halfZ);

    /// Forward intersection test
    Vector3 before = transform * Vector3(-50_cm, -1_m, -1_m);
    Vector3 onit1 = transform * Vector3(0_m, 0_m, 0_m);
    Vector3 onitP = transform * Vector3(1_cm, 0_m, 23_um);
    Vector3 after = transform * Vector3(33_cm, 12_mm, 1_m);
    Vector3 outside = transform * Vector3(2., 0., 100_m);

    // Linear transform
    auto lTransform = transform.linear();
    Vector3 direction = lTransform * Vector3(2_cm, 3_cm, 5_cm).normalized();
    Vector3 normalP = lTransform * Vector3(0, 1., 0.).normalized();
    Vector3 parallel = lTransform * Vector3(0, 0., 1.).normalized();

    // A random intersection form backward
    // Intersect forward
    auto fIntersection = aLine->intersect(tgContext, before, direction, true);
    // The intersection MUST be valid
    BOOST_CHECK(fIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(fIntersection[0].status() == Intersection3D::Status::reachable);
    // The path length MUST be positive
    BOOST_CHECK(fIntersection[0].pathLength() > 0.);
    // The intersection MUST be unique
    BOOST_CHECK(!fIntersection[1]);

    // On surface intersection - on the straw with random direction
    auto oIntersection = aLine->intersect(tgContext, onit1, direction, true);
    // The intersection MUST be valid
    BOOST_CHECK(oIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(oIntersection[0].status() == Intersection3D::Status::onSurface);
    // The path length MUST be positive
    BOOST_CHECK(std::abs(oIntersection[0].pathLength()) < s_onSurfaceTolerance);
    // The intersection MUST be unique
    BOOST_CHECK(!oIntersection[1]);

    // On surface intersecion - on the surface with normal vector
    oIntersection = aLine->intersect(tgContext, onitP, normalP, true);
    // The intersection MUST be valid
    BOOST_CHECK(oIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(oIntersection[0].status() == Intersection3D::Status::onSurface);
    // The path length MUST be positive
    BOOST_CHECK(std::abs(oIntersection[0].pathLength()) < s_onSurfaceTolerance);
    // The intersection MUST be unique
    BOOST_CHECK(!oIntersection[1]);

    // Intersect backwards
    auto bIntersection = aLine->intersect(tgContext, after, direction, true);
    // The intersection MUST be valid
    BOOST_CHECK(bIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(bIntersection[0].status() == Intersection3D::Status::reachable);
    // The path length MUST be negative
    BOOST_CHECK(bIntersection[0].pathLength() < 0.);
    // The intersection MUST be unique
    BOOST_CHECK(!bIntersection[1]);

    // An out of bounds attempt: missed
    auto mIntersection = aLine->intersect(tgContext, outside, direction, true);
    // The intersection MUST NOT be valid
    BOOST_CHECK(!mIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(mIntersection[0].status() == Intersection3D::Status::missed);
    // The path length MUST be negative
    BOOST_CHECK(mIntersection[0].pathLength() < 0.);
    // The intersection MUST be unique
    BOOST_CHECK(!mIntersection[1]);

    // An invalid attempt
    auto iIntersection = aLine->intersect(tgContext, before, parallel, true);
    // The intersection MUST NOT be valid
    BOOST_CHECK(!iIntersection[0]);
    // The intersection MUST be reachable
    BOOST_CHECK(iIntersection[0].status() ==
                Intersection3D::Status::unreachable);
    // The intersection MUST be unique
    BOOST_CHECK(!iIntersection[1]);
  };

  // In a nominal world
  testLineAppraoch(Transform3::Identity());

  // In a system somewhere away
  testLineAppraoch(aTransform);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
