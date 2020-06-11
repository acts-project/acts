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

BOOST_AUTO_TEST_CASE(PlaneSurfaceEqualityOperators) {
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

/// Unit test for testing PlaneSurface extent via Polyhedron representation
BOOST_AUTO_TEST_CASE(PlaneSurfaceExtent) {
  // First test - non-rotated
  static const Transform3D planeZX =
      AngleAxis3D(-0.5 * M_PI, Vector3D::UnitX()) *
      AngleAxis3D(-0.5 * M_PI, Vector3D::UnitZ()) * Transform3D::Identity();

  double rHx = 2.;
  double rHy = 4.;
  double yPs = 3.;
  auto rBounds = std::make_shared<RectangleBounds>(rHx, rHy);

  auto plane = Surface::makeShared<PlaneSurface>(
      std::make_shared<Transform3D>(Translation3D(Vector3D(0., yPs, 0.)) *
                                    planeZX),
      rBounds);

  auto planeExtent = plane->polyhedronRepresentation(tgContext, 1).extent();

  CHECK_CLOSE_ABS(planeExtent.min(binZ), -rHx, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.max(binZ), rHx, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.min(binX), -rHy, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.max(binX), rHy, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.min(binY), yPs, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.max(binY), yPs, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.min(binR), yPs, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtent.max(binR), std::sqrt(yPs * yPs + rHy * rHy),
                  s_onSurfaceTolerance);

  // Now rotate
  double alpha = 0.123;
  auto planeRot = Surface::makeShared<PlaneSurface>(
      std::make_shared<Transform3D>(Translation3D(Vector3D(0., yPs, 0.)) *
                                    AngleAxis3D(alpha, Vector3D(0., 0., 1.)) *
                                    planeZX),
      rBounds);

  auto planeExtentRot =
      planeRot->polyhedronRepresentation(tgContext, 1).extent();
  CHECK_CLOSE_ABS(planeExtentRot.min(binZ), -rHx, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.max(binZ), rHx, s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.min(binX), -rHy * std::cos(alpha),
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.max(binX), rHy * std::cos(alpha),
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.min(binY), yPs - rHy * std::sin(alpha),
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.max(binY), yPs + rHy * std::sin(alpha),
                  s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(planeExtentRot.min(binR), yPs * std::cos(alpha),
                  s_onSurfaceTolerance);
}

/// Unit test for testing PlaneSurface alignment derivatives
BOOST_AUTO_TEST_CASE(PlaneSurfaceAlignment) {
  // bounds object, rectangle type
  auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);
  /// Test clone method
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto planeSurfaceObject =
      Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  const auto& rotation = pTransform->rotation();
  // The local frame z axis
  const Vector3D localZAxis = rotation.col(2);
  // Check the local z axis is aligned to global z axis
  CHECK_CLOSE_ABS(localZAxis, Vector3D(0., 0., 1.), 1e-15);

  /// Define the track (local) position and direction
  Vector2D localPosition{1, 2};
  Vector3D momentum{0, 0, 1};
  Vector3D direction = momentum.normalized();
  /// Get the global position
  Vector3D globalPosition{0, 0, 0};
  planeSurfaceObject->localToGlobal(tgContext, localPosition, momentum,
                                    globalPosition);

  // Call the function to calculate the derivative of local frame axes w.r.t its
  // rotation
  const auto& [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);

  // (a) Test the derivative of path length w.r.t. alignment parameters
  const AlignmentRowVector& alignToPath =
      planeSurfaceObject->alignmentToPathDerivative(tgContext, rotToLocalZAxis,
                                                    globalPosition, direction);
  // The expected results
  AlignmentRowVector expAlignToPath = AlignmentRowVector::Zero();
  expAlignToPath << 0, 0, 1, 2, -1, 0;
  // Check if the calculated derivative is as expected
  CHECK_CLOSE_ABS(alignToPath, expAlignToPath, 1e-10);

  // (b) Test the derivative of bound track parameters local position w.r.t.
  // position in local 3D Cartesian coordinates
  const auto& loc3DToLocBound =
      planeSurfaceObject->localCartesianToBoundLocalDerivative(tgContext,
                                                               globalPosition);
  // For plane surface, this should be identity matrix
  CHECK_CLOSE_ABS(loc3DToLocBound, LocalCartesianToBoundLocalMatrix::Identity(),
                  1e-10);

  // (c) Test the derivative of bound parameters (only test loc0, loc1 here)
  // w.r.t. alignment parameters
  FreeVector derivatives = FreeVector::Zero();
  derivatives.segment<3>(0) = momentum;
  const AlignmentToBoundMatrix& alignToBound =
      planeSurfaceObject->alignmentToBoundDerivative(tgContext, derivatives,
                                                     globalPosition, direction);
  const AlignmentRowVector& alignToloc0 = alignToBound.block<1, 6>(0, 0);
  const AlignmentRowVector& alignToloc1 = alignToBound.block<1, 6>(1, 0);
  // The expected results
  AlignmentRowVector expAlignToloc0;
  expAlignToloc0 << -1, 0, 0, 0, 0, 2;
  AlignmentRowVector expAlignToloc1;
  expAlignToloc1 << 0, -1, 0, 0, 0, -1;
  // Check if the calculated derivatives are as expected
  CHECK_CLOSE_ABS(alignToloc0, expAlignToloc0, 1e-10);
  CHECK_CLOSE_ABS(alignToloc1, expAlignToloc1, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
