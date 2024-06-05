// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <memory>
#include <ostream>
#include <string>
#include <utility>

namespace Acts {
class AssertionFailureException;
}  // namespace Acts

namespace Acts::Test {

// Create a test context
GeometryContext testContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(CylinderSurfaces)
/// Unit test for creating compliant/non-compliant CylinderSurface object
BOOST_AUTO_TEST_CASE(CylinderSurfaceConstruction) {
  // CylinderSurface default constructor is deleted
  //
  /// Constructor with transform, radius and halfZ
  double radius(1.0), halfZ(10.), halfPhiSector(M_PI / 8.);
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
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
      testContext, *cylinderSurfaceObject, pTransform);
  BOOST_CHECK_EQUAL(copiedTransformedCylinderSurface->type(),
                    Surface::Cylinder);

  /// Construct with nullptr bounds
  BOOST_CHECK_THROW(auto nullBounds = Surface::makeShared<CylinderSurface>(
                        Transform3::Identity(), nullptr),
                    AssertionFailureException);
}
//
/// Unit test for testing CylinderSurface properties
BOOST_AUTO_TEST_CASE(CylinderSurfaceProperties) {
  /// Test clone method
  double radius(1.0), halfZ(10.);
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  auto cylinderSurfaceObject =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
  //
  /// Test type (redundant)
  BOOST_CHECK_EQUAL(cylinderSurfaceObject->type(), Surface::Cylinder);
  //
  /// Test binningPosition
  Vector3 binningPosition{0., 1., 2.};
  CHECK_CLOSE_ABS(
      cylinderSurfaceObject->binningPosition(testContext, BinningValue::binPhi),
      binningPosition, 1e-9);
  //
  /// Test referenceFrame
  double rootHalf = std::sqrt(0.5);
  Vector3 globalPosition{rootHalf, 1. - rootHalf, 0.};
  Vector3 globalPositionZ{rootHalf, 1. - rootHalf, 2.0};
  Vector3 momentum{15., 15., 15.};
  Vector3 momentum2{6.6, -3., 2.};
  RotationMatrix3 expectedFrame;
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
  Vector3 origin{0., 0., 0.};
  Vector3 normal3D = {0., -1., 0.};
  CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(testContext, origin), normal3D,
                  1e-9);

  Vector3 pos45deg = {rootHalf, 1 + rootHalf, 0.};
  Vector3 pos45degZ = {rootHalf, 1 + rootHalf, 4.};
  Vector3 normal45deg = {rootHalf, rootHalf, 0.};
  // test the normal vector
  CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(testContext, pos45deg),
                  normal45deg, 1e-6 * rootHalf);
  // test that the normal vector is independent of z coordinate
  CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(testContext, pos45degZ),
                  normal45deg, 1e-6 * rootHalf);
  //
  /// Test normal given 2D rphi position
  Vector2 positionPiBy2(1.0, 0.);
  Vector3 normalAtPiBy2{std::cos(1.), std::sin(1.), 0.};
  CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(testContext, positionPiBy2),
                  normalAtPiBy2, 1e-9);

  //
  /// Test rotational symmetry axis
  Vector3 symmetryAxis{0., 0., 1.};
  CHECK_CLOSE_ABS(cylinderSurfaceObject->rotSymmetryAxis(testContext),
                  symmetryAxis, 1e-9);
  //
  /// Test bounds
  BOOST_CHECK_EQUAL(cylinderSurfaceObject->bounds().type(),
                    SurfaceBounds::eCylinder);
  //
  /// Test localToGlobal
  Vector2 localPosition{0., 0.};
  globalPosition = cylinderSurfaceObject->localToGlobal(
      testContext, localPosition, momentum);
  Vector3 expectedPosition{1, 1, 2};
  BOOST_CHECK_EQUAL(globalPosition, expectedPosition);
  //
  /// Testing globalToLocal
  localPosition = cylinderSurfaceObject
                      ->globalToLocal(testContext, globalPosition, momentum)
                      .value();
  Vector2 expectedLocalPosition{0., 0.};
  BOOST_CHECK_EQUAL(localPosition, expectedLocalPosition);
  //
  /// Test isOnSurface
  Vector3 offSurface{100, 1, 2};
  BOOST_CHECK(cylinderSurfaceObject->isOnSurface(
      testContext, globalPosition, momentum, BoundaryCheck(true)));
  BOOST_CHECK(!cylinderSurfaceObject->isOnSurface(
      testContext, offSurface, momentum, BoundaryCheck(true)));
  //
  /// intersection test
  Vector3 direction{-1., 0, 0};
  auto sfIntersection = cylinderSurfaceObject->intersect(
      testContext, offSurface, direction, BoundaryCheck(false));
  Intersection3D expectedIntersect{Vector3{1, 1, 2}, 99.,
                                   Intersection3D::Status::reachable};
  BOOST_CHECK(sfIntersection[0]);
  CHECK_CLOSE_ABS(sfIntersection[0].position(), expectedIntersect.position(),
                  1e-9);
  CHECK_CLOSE_ABS(sfIntersection[0].pathLength(),
                  expectedIntersect.pathLength(), 1e-9);
  // there is a second solution & and it should be valid
  BOOST_CHECK(sfIntersection[1]);
  // And it's path should be further away then the primary solution
  double pn = sfIntersection[0].pathLength();
  double pa = sfIntersection[1].pathLength();
  BOOST_CHECK_LT(std::abs(pn), std::abs(pa));
  BOOST_CHECK_EQUAL(sfIntersection.object(), cylinderSurfaceObject.get());

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
     Bounds  : Acts::CylinderBounds: (radius, halfLengthZ, halfPhiSector, averagePhi, bevelMinZ, bevelMaxZ) = (1.0000000, 10.0000000, 3.1415927, 0.0000000, 0.0000000, 0.0000000)"));
}

BOOST_AUTO_TEST_CASE(CylinderSurfaceEqualityOperators) {
  double radius(1.0), halfZ(10.);
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
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
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 6.6, 5.4);
  *assignedCylinderSurface = *cylinderSurfaceObject;
  /// Test equality of assigned to original
  BOOST_CHECK(*assignedCylinderSurface == *cylinderSurfaceObject);
}

/// Unit test for testing CylinderSurface properties
BOOST_AUTO_TEST_CASE(CylinderSurfaceExtent) {
  // Some radius and half length
  double radius(1.0), halfZ(10.);
  Translation3 translation{0., 0., 2.};
  auto pTransform = Transform3(translation);
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
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  auto cylinderSurfaceObject =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);

  const auto& rotation = pTransform.rotation();
  // The local frame z axis
  const Vector3 localZAxis = rotation.col(2);
  // Check the local z axis is aligned to global z axis
  CHECK_CLOSE_ABS(localZAxis, Vector3(0., 0., 1.), 1e-15);

  /// Define the track (global) position and direction
  Vector3 globalPosition{0, 2, 2};

  // Test the derivative of bound track parameters local position w.r.t.
  // position in local 3D Cartesian coordinates
  const auto& loc3DToLocBound =
      cylinderSurfaceObject->localCartesianToBoundLocalDerivative(
          testContext, globalPosition);
  // Check if the result is as expected
  ActsMatrix<2, 3> expLoc3DToLocBound = ActsMatrix<2, 3>::Zero();
  expLoc3DToLocBound << -1, 0, 0, 0, 0, 1;
  CHECK_CLOSE_ABS(loc3DToLocBound, expLoc3DToLocBound, 1e-10);
}

BOOST_AUTO_TEST_CASE(CylinderSurfaceBinningPosition) {
  using namespace Acts::UnitLiterals;
  Vector3 s{5_mm, 7_mm, 10_cm};
  Transform3 trf;
  trf = Translation3(s) * AngleAxis3{0.5, Vector3::UnitZ()};

  double r = 300;
  double halfZ = 330;
  double averagePhi = 0.1;

  auto bounds =
      std::make_shared<CylinderBounds>(r, halfZ, M_PI / 8, averagePhi);
  auto cylinder = Acts::Surface::makeShared<CylinderSurface>(trf, bounds);

  Vector3 exp = Vector3{r * std::cos(averagePhi), r * std::sin(averagePhi), 0};
  exp = trf * exp;

  Vector3 bp = cylinder->binningPosition(testContext, binR);
  CHECK_CLOSE_ABS(bp, exp, 1e-10);
  CHECK_CLOSE_ABS(cylinder->binningPositionValue(testContext, binR),
                  VectorHelpers::perp(exp), 1e-10);

  bp = cylinder->binningPosition(testContext, binRPhi);
  CHECK_CLOSE_ABS(bp, exp, 1e-10);
  CHECK_CLOSE_ABS(cylinder->binningPositionValue(testContext, binRPhi),
                  VectorHelpers::phi(exp) * VectorHelpers::perp(exp), 1e-10);

  for (auto b : {binX, binY, binZ, binEta, binH, binMag}) {
    BOOST_TEST_CONTEXT("binValue: " << b) {
      BOOST_CHECK_EQUAL(cylinder->binningPosition(testContext, b),
                        cylinder->center(testContext));
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
