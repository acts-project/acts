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

#include <limits>

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace tt = boost::test_tools;
using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(PlaneSurfaces)
/// Unit test for creating compliant/non-compliant PlaneSurface object
BOOST_AUTO_TEST_CASE(PlaneSurfaceConstruction) {
  // PlaneSurface default constructor is deleted
  // bounds object, rectangle type
  auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);
  /// Constructor with transform pointer, null or valid, alpha and symmetry
  /// indicator
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  std::shared_ptr<const Transform3D> pNullTransform{};
  // constructor with nullptr transform
  BOOST_CHECK_EQUAL(
      Surface::makeShared<PlaneSurface>(pNullTransform, rBounds)->type(),
      Surface::Plane);
  // constructor with transform
  BOOST_CHECK_EQUAL(
      Surface::makeShared<PlaneSurface>(pTransform, rBounds)->type(),
      Surface::Plane);
  /// Copy constructor
  auto planeSurfaceObject =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  auto copiedPlaneSurface =
      Surface::makeShared<PlaneSurface>(*planeSurfaceObject);
  BOOST_CHECK_EQUAL(copiedPlaneSurface->type(), Surface::Plane);
  BOOST_CHECK(*copiedPlaneSurface == *planeSurfaceObject);
  //
  /// Copied and transformed
  auto copiedTransformedPlaneSurface = Surface::makeShared<PlaneSurface>(
      tgContext, *planeSurfaceObject, *pTransform);
  BOOST_CHECK_EQUAL(copiedTransformedPlaneSurface->type(), Surface::Plane);

  /// Construct with nullptr bounds
  DetectorElementStub detElem;
  BOOST_CHECK_THROW(
      auto nullBounds = Surface::makeShared<PlaneSurface>(nullptr, detElem),
      AssertionFailureException);
}
//
/// Unit test for testing PlaneSurface properties
BOOST_AUTO_TEST_CASE(PlaneSurfaceProperties) {
  // bounds object, rectangle type
  auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);
  /// Test clone method
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto planeSurfaceObject =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  // Is it in the right place?
  Translation3D translation2{0., 2., 4.};
  auto pTransform2 = std::make_shared<const Transform3D>(translation2);
  auto planeSurfaceObject2 =
      Surface::makeShared<PlaneSurface>(pTransform2, rBounds);
  /// Test type (redundant)
  BOOST_CHECK_EQUAL(planeSurfaceObject->type(), Surface::Plane);
  //
  /// Test binningPosition
  Vector3D binningPosition{0., 1., 2.};
  BOOST_CHECK_EQUAL(
      planeSurfaceObject->binningPosition(tgContext, BinningValue::binX),
      binningPosition);
  //
  /// Test referenceFrame
  Vector3D globalPosition{2.0, 2.0, 0.0};
  Vector3D momentum{1.e6, 1.e6, 1.e6};
  RotationMatrix3D expectedFrame;
  expectedFrame << 1., 0., 0., 0., 1., 0., 0., 0., 1.;

  CHECK_CLOSE_OR_SMALL(
      planeSurfaceObject->referenceFrame(tgContext, globalPosition, momentum),
      expectedFrame, 1e-6, 1e-9);
  //
  /// Test normal, given 3D position
  Vector3D normal3D(0., 0., 1.);
  BOOST_CHECK_EQUAL(planeSurfaceObject->normal(tgContext), normal3D);
  //
  /// Test bounds
  BOOST_CHECK_EQUAL(planeSurfaceObject->bounds().type(),
                    SurfaceBounds::eRectangle);

  /// Test localToGlobal
  Vector2D localPosition{1.5, 1.7};
  planeSurfaceObject->localToGlobal(tgContext, localPosition, momentum,
                                    globalPosition);
  //
  // expected position is the translated one
  Vector3D expectedPosition{1.5 + translation.x(), 1.7 + translation.y(),
                            translation.z()};

  CHECK_CLOSE_REL(globalPosition, expectedPosition, 1e-2);
  //
  /// Testing globalToLocal
  planeSurfaceObject->globalToLocal(tgContext, globalPosition, momentum,
                                    localPosition);
  Vector2D expectedLocalPosition{1.5, 1.7};

  CHECK_CLOSE_REL(localPosition, expectedLocalPosition, 1e-2);

  /// Test isOnSurface
  Vector3D offSurface{0, 1, -2.};
  BOOST_CHECK(planeSurfaceObject->isOnSurface(tgContext, globalPosition,
                                              momentum, true));
  BOOST_CHECK(
      !planeSurfaceObject->isOnSurface(tgContext, offSurface, momentum, true));
  //
  /// intersectionEstimate
  Vector3D direction{0., 0., 1.};
  auto intersect = planeSurfaceObject->intersectionEstimate(
      tgContext, offSurface, direction, true);
  Intersection expectedIntersect{Vector3D{0, 1, 2}, 4.,
                                 Intersection::Status::reachable};
  BOOST_CHECK(bool(intersect));
  BOOST_CHECK_EQUAL(intersect.position, expectedIntersect.position);
  BOOST_CHECK_EQUAL(intersect.pathLength, expectedIntersect.pathLength);
  //

  /// Test pathCorrection
  CHECK_CLOSE_REL(planeSurfaceObject->pathCorrection(tgContext, offSurface,
                                                     momentum.normalized()),
                  std::sqrt(3), 0.01);
  //
  /// Test name
  BOOST_CHECK_EQUAL(planeSurfaceObject->name(),
                    std::string("Acts::PlaneSurface"));
  //
  /// Test dump
  // TODO 2017-04-12 msmk: check how to correctly check output
  //    boost::test_tools::output_test_stream dumpOuput;
  //    planeSurfaceObject.toStream(dumpOuput);
  //    BOOST_CHECK(dumpOuput.is_equal(
  //      "Acts::PlaneSurface\n"
  //      "    Center position  (x, y, z) = (0.0000, 1.0000, 2.0000)\n"
  //      "    Rotation:             colX = (1.000000, 0.000000, 0.000000)\n"
  //      "                          colY = (0.000000, 1.000000, 0.000000)\n"
  //      "                          colZ = (0.000000, 0.000000, 1.000000)\n"
  //      "    Bounds  : Acts::ConeBounds: (tanAlpha, minZ, maxZ, averagePhi,
  //      halfPhiSector) = (0.4142136, 0.0000000, inf, 0.0000000,
  //      3.1415927)"));
}

BOOST_AUTO_TEST_CASE(EqualityOperators) {
  // rectangle bounds
  auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto planeSurfaceObject =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  auto planeSurfaceObject2 =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  //
  /// Test equality operator
  BOOST_CHECK(*planeSurfaceObject == *planeSurfaceObject2);
  //
  BOOST_TEST_CHECKPOINT(
      "Create and then assign a PlaneSurface object to the existing one");
  /// Test assignment
  auto assignedPlaneSurface =
      Surface::makeShared<PlaneSurface>(nullptr, nullptr);
  *assignedPlaneSurface = *planeSurfaceObject;
  /// Test equality of assigned to original
  BOOST_CHECK(*assignedPlaneSurface == *planeSurfaceObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
