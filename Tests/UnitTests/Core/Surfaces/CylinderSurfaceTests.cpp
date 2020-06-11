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

#include "Acts/Surfaces/CylinderSurface.hpp"
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
GeometryContext testContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(CylinderSurfaces)
/// Unit test for creating compliant/non-compliant CylinderSurface object
BOOST_AUTO_TEST_CASE(CylinderSurfaceConstruction) {
  // CylinderSurface default constructor is deleted
  //
  /// Constructor with transform pointer, null or valid, radius and halfZ
  double radius(1.0), halfZ(10.), halfPhiSector(M_PI / 8.);
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  std::shared_ptr<const Transform3D> pNullTransform{};
  BOOST_CHECK_EQUAL(
      Surface::makeShared<CylinderSurface>(pNullTransform, radius, halfZ)
          ->type(),
      Surface::Cylinder);
  BOOST_CHECK_EQUAL(
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ)->type(),
      Surface::Cylinder);
  //
  /// Constructor with transform pointer, radius, halfZ and halfPhiSector
  BOOST_CHECK_EQUAL(Surface::makeShared<CylinderSurface>(pTransform, radius,
                                                         halfZ, halfPhiSector)
                        ->type(),
                    Surface::Cylinder);

  /// Constructor with transform and CylinderBounds pointer
  auto pCylinderBounds = std::make_shared<const CylinderBounds>(radius, halfZ);
  BOOST_CHECK_EQUAL(
      Surface::makeShared<CylinderSurface>(pTransform, pCylinderBounds)->type(),
      Surface::Cylinder);
  //
  //
  /// Copy constructor
  auto cylinderSurfaceObject =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
  auto copiedCylinderSurface =
      Surface::makeShared<CylinderSurface>(*cylinderSurfaceObject);
  BOOST_CHECK_EQUAL(copiedCylinderSurface->type(), Surface::Cylinder);
  BOOST_CHECK(*copiedCylinderSurface == *cylinderSurfaceObject);
  //
  /// Copied and transformed
  auto copiedTransformedCylinderSurface = Surface::makeShared<CylinderSurface>(
      testContext, *cylinderSurfaceObject, *pTransform);
  BOOST_CHECK_EQUAL(copiedTransformedCylinderSurface->type(),
                    Surface::Cylinder);

  /// Construct with nullptr bounds
  BOOST_CHECK_THROW(
      auto nullBounds = Surface::makeShared<CylinderSurface>(nullptr, nullptr),
      AssertionFailureException);
}
//
/// Unit test for testing CylinderSurface properties
BOOST_AUTO_TEST_CASE(CylinderSurfaceProperties) {
  /// Test clone method
  double radius(1.0), halfZ(10.);
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto cylinderSurfaceObject =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
  //
  /// Test type (redundant)
  BOOST_CHECK_EQUAL(cylinderSurfaceObject->type(), Surface::Cylinder);
  //
  /// Test binningPosition
  Vector3D binningPosition{0., 1., 2.};
  CHECK_CLOSE_ABS(
      cylinderSurfaceObject->binningPosition(testContext, BinningValue::binPhi),
      binningPosition, 1e-9);
  //
  /// Test referenceFrame
  double rootHalf = std::sqrt(0.5);
  Vector3D globalPosition{rootHalf, 1. - rootHalf, 0.};
  Vector3D globalPositionZ{rootHalf, 1. - rootHalf, 2.0};
  Vector3D momentum{15., 15., 15.};
  Vector3D momentum2{6.6, -3., 2.};
  RotationMatrix3D expectedFrame;
  expectedFrame << rootHalf, 0., rootHalf, rootHalf, 0., -rootHalf, 0., 1., 0.;
  // check without shift
  CHECK_CLOSE_OR_SMALL(cylinderSurfaceObject->referenceFrame(
                           testContext, globalPosition, momentum),
                       expectedFrame, 1e-6, 1e-9);
  // check with shift and different momentum
  CHECK_CLOSE_OR_SMALL(cylinderSurfaceObject->referenceFrame(
                           testContext, globalPositionZ, momentum2),
                       expectedFrame, 1e-6, 1e-9);
  //
  /// Test normal, given 3D position
  Vector3D origin{0., 0., 0.};
  Vector3D normal3D = {0., -1., 0.};
  CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(testContext, origin), normal3D,
                  1e-9);

  Vector3D pos45deg = {rootHalf, 1 + rootHalf, 0.};
  Vector3D pos45degZ = {rootHalf, 1 + rootHalf, 4.};
  Vector3D normal45deg = {rootHalf, rootHalf, 0.};
  // test the normal vector
  CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(testContext, pos45deg),
                  normal45deg, 1e-6 * rootHalf);
  // thest that the normal vector is independent of z coordinate
  CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(testContext, pos45degZ),
                  normal45deg, 1e-6 * rootHalf);
  //
  /// Test normal given 2D rphi position
  Vector2D positionPiBy2(1.0, 0.);
  Vector3D normalAtPiBy2{std::cos(1.), std::sin(1.), 0.};
  CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(testContext, positionPiBy2),
                  normalAtPiBy2, 1e-9);

  //
  /// Test rotational symmetry axis
  Vector3D symmetryAxis{0., 0., 1.};
  CHECK_CLOSE_ABS(cylinderSurfaceObject->rotSymmetryAxis(testContext),
                  symmetryAxis, 1e-9);
  //
  /// Test bounds
  BOOST_CHECK_EQUAL(cylinderSurfaceObject->bounds().type(),
                    SurfaceBounds::eCylinder);
  //
  /// Test localToGlobal
  Vector2D localPosition{0., 0.};
  cylinderSurfaceObject->localToGlobal(testContext, localPosition, momentum,
                                       globalPosition);
  Vector3D expectedPosition{1, 1, 2};
  BOOST_CHECK_EQUAL(globalPosition, expectedPosition);
  //
  /// Testing globalToLocal
  cylinderSurfaceObject->globalToLocal(testContext, globalPosition, momentum,
                                       localPosition);
  Vector2D expectedLocalPosition{0., 0.};
  BOOST_CHECK_EQUAL(localPosition, expectedLocalPosition);
  //
  /// Test isOnSurface
  Vector3D offSurface{100, 1, 2};
  BOOST_CHECK(cylinderSurfaceObject->isOnSurface(testContext, globalPosition,
                                                 momentum, true));
  BOOST_CHECK(!cylinderSurfaceObject->isOnSurface(testContext, offSurface,
                                                  momentum, true));
  //
  /// intersectionEstimate
  Vector3D direction{-1., 0, 0};
  auto intersect = cylinderSurfaceObject->intersectionEstimate(
      testContext, offSurface, direction, false);
  Intersection expectedIntersect{Vector3D{1, 1, 2}, 99.,
                                 Intersection::Status::reachable};
  // check the result
  BOOST_CHECK(bool(intersect));
  CHECK_CLOSE_ABS(intersect.position, expectedIntersect.position, 1e-9);
  CHECK_CLOSE_ABS(intersect.pathLength, expectedIntersect.pathLength, 1e-9);

  /// intersect
  auto surfaceIntersect = cylinderSurfaceObject->intersect(
      testContext, offSurface, direction, false);
  BOOST_CHECK(bool(surfaceIntersect));
  // there is a second solution & and it should be valid
  BOOST_CHECK(surfaceIntersect.alternative);
  // And it's path should be further away then the primary solution
  double pn = surfaceIntersect.intersection.pathLength;
  double pa = surfaceIntersect.alternative.pathLength;
  BOOST_CHECK(pn * pn < pa * pa);
  //
  /// Test pathCorrection
  CHECK_CLOSE_REL(cylinderSurfaceObject->pathCorrection(testContext, offSurface,
                                                        momentum.normalized()),
                  std::sqrt(3.), 0.01);
  //
  /// Test name
  BOOST_CHECK_EQUAL(cylinderSurfaceObject->name(),
                    std::string("Acts::CylinderSurface"));
  //
  /// Test dump
  boost::test_tools::output_test_stream dumpOuput;
  cylinderSurfaceObject->toStream(testContext, dumpOuput);
  BOOST_CHECK(
      dumpOuput.is_equal("Acts::CylinderSurface\n\
     Center position  (x, y, z) = (0.0000, 1.0000, 2.0000)\n\
     Rotation:             colX = (1.000000, 0.000000, 0.000000)\n\
                           colY = (0.000000, 1.000000, 0.000000)\n\
                           colZ = (0.000000, 0.000000, 1.000000)\n\
     Bounds  : Acts::CylinderBounds: (radius, halfLengthZ, halfPhiSector, averagePhi) = (1.0000000, 10.0000000, 3.1415927, 0.0000000)"));
}

BOOST_AUTO_TEST_CASE(CylinderSurfaceEqualityOperators) {
  double radius(1.0), halfZ(10.);
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto cylinderSurfaceObject =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
  //
  auto cylinderSurfaceObject2 =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
  //
  /// Test equality operator
  BOOST_CHECK(*cylinderSurfaceObject == *cylinderSurfaceObject2);
  //
  BOOST_TEST_CHECKPOINT(
      "Create and then assign a CylinderSurface object to the existing one");
  /// Test assignment
  auto assignedCylinderSurface =
      Surface::makeShared<CylinderSurface>(nullptr, 6.6, 5.4);
  *assignedCylinderSurface = *cylinderSurfaceObject;
  /// Test equality of assigned to original
  BOOST_CHECK(*assignedCylinderSurface == *cylinderSurfaceObject);
}

/// Unit test for testing CylinderSurface properties
BOOST_AUTO_TEST_CASE(CylinderSurfaceExtent) {
  // Some radius and half length
  double radius(1.0), halfZ(10.);
  Translation3D translation{0., 0., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto cylinderSurface =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
  // The Extent, let's measure it
  auto cylinderExtent =
      cylinderSurface->polyhedronRepresentation(testContext, 1).extent();

  CHECK_CLOSE_ABS(-8, cylinderExtent.min(binZ), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(12, cylinderExtent.max(binZ), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(radius, cylinderExtent.min(binR), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(radius, cylinderExtent.max(binR), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(-radius, cylinderExtent.min(binX), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(radius, cylinderExtent.max(binX), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(-radius, cylinderExtent.min(binY), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(radius, cylinderExtent.max(binY), s_onSurfaceTolerance);
}

/// Unit test for testing CylinderSurface alignment derivatives
BOOST_AUTO_TEST_CASE(CylinderSurfaceAlignment) {
  double radius(1.0), halfZ(10.);
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto cylinderSurfaceObject =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);

  const auto& rotation = pTransform->rotation();
  // The local frame z axis
  const Vector3D localZAxis = rotation.col(2);
  // Check the local z axis is aligned to global z axis
  CHECK_CLOSE_ABS(localZAxis, Vector3D(0., 0., 1.), 1e-15);

  /// Define the track (global) position and direction
  Vector3D globalPosition{0, 2, 2};

  // Test the derivative of bound track parameters local position w.r.t.
  // position in local 3D Cartesian coordinates
  const auto& loc3DToLocBound =
      cylinderSurfaceObject->localCartesianToBoundLocalDerivative(
          testContext, globalPosition);
  // Check if the result is as expected
  LocalCartesianToBoundLocalMatrix expLoc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Zero();
  expLoc3DToLocBound << -1, 0, 0, 0, 0, 1;
  CHECK_CLOSE_ABS(loc3DToLocBound, expLoc3DToLocBound, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
