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

#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

namespace Acts {

namespace Test {
// using boost::test_tools::output_test_stream;
// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit tests for creating DiscSurface object
BOOST_AUTO_TEST_CASE(DiscSurfaceConstruction) {
  // default constructor is deleted
  // scaffolding...
  double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
  //
  /// Test DiscSurface fully specified constructor but no transform
  BOOST_CHECK_NO_THROW(
      Surface::makeShared<DiscSurface>(nullptr, rMin, rMax, halfPhiSector));
  //
  /// Test DiscSurface constructor with default halfPhiSector
  BOOST_CHECK_NO_THROW(Surface::makeShared<DiscSurface>(nullptr, rMin, rMax));
  //
  /// Test DiscSurface constructor with a transform specified
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  BOOST_CHECK_NO_THROW(
      Surface::makeShared<DiscSurface>(pTransform, rMin, rMax, halfPhiSector));
  //
  /// Copy constructed DiscSurface
  auto anotherDiscSurface =
      Surface::makeShared<DiscSurface>(pTransform, rMin, rMax, halfPhiSector);
  // N.B. Just using
  // BOOST_CHECK_NO_THROW(Surface::makeShared<DiscSurface>(anotherDiscSurface))
  // tries to call
  // the (deleted) default constructor.
  auto copiedSurface = Surface::makeShared<DiscSurface>(*anotherDiscSurface);
  BOOST_TEST_MESSAGE("Copy constructed DiscSurface ok");
  //
  /// Copied and transformed DiscSurface
  BOOST_CHECK_NO_THROW(Surface::makeShared<DiscSurface>(
      tgContext, *anotherDiscSurface, *pTransform));

  /// Construct with nullptr bounds
  DetectorElementStub detElem;
  BOOST_CHECK_THROW(
      auto nullBounds = Surface::makeShared<DiscSurface>(nullptr, detElem),
      AssertionFailureException);
}

/// Unit tests of all named methods
BOOST_AUTO_TEST_CASE(DiscSurfaceProperties, *utf::expected_failures(2)) {
  Vector3D origin3D{0, 0, 0};
  std::shared_ptr<const Transform3D> pTransform;  // nullptr
  double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
  auto discSurfaceObject =
      Surface::makeShared<DiscSurface>(pTransform, rMin, rMax, halfPhiSector);
  //
  /// Test type
  BOOST_CHECK_EQUAL(discSurfaceObject->type(), Surface::Disc);
  //
  /// Test normal, no local position specified
  Vector3D zAxis{0, 0, 1};
  BOOST_CHECK_EQUAL(discSurfaceObject->normal(tgContext), zAxis);
  //
  /// Test normal, local position specified
  Vector2D lpos(2.0, 0.05);
  BOOST_CHECK_EQUAL(discSurfaceObject->normal(tgContext, lpos), zAxis);
  //
  /// Test binningPosition
  // auto binningPosition=
  // discSurfaceObject.binningPosition(BinningValue::binRPhi );
  // std::cout<<binningPosition<<std::endl;
  BOOST_CHECK_EQUAL(
      discSurfaceObject->binningPosition(tgContext, BinningValue::binRPhi),
      origin3D);
  //
  /// Test bounds
  BOOST_CHECK_EQUAL(discSurfaceObject->bounds().type(), SurfaceBounds::eDisc);
  //
  Vector3D ignoredMomentum{0., 0., 0.};
  /// Test isOnSurface()
  Vector3D point3DNotInSector{0.0, 1.2, 0};
  Vector3D point3DOnSurface{1.2, 0.0, 0};
  BOOST_CHECK(!discSurfaceObject->isOnSurface(
      tgContext, point3DNotInSector, ignoredMomentum, true));  // passes
  BOOST_CHECK(discSurfaceObject->isOnSurface(tgContext, point3DOnSurface,
                                             ignoredMomentum, true));  // passes
  //
  /// Test localToGlobal
  Vector3D returnedPosition{10.9, 8.7, 6.5};
  Vector3D expectedPosition{1.2, 0, 0};
  Vector2D rPhiOnDisc{1.2, 0.0};
  Vector2D rPhiNotInSector{1.2, M_PI};  // outside sector at Phi=0, +/- pi/8
  discSurfaceObject->localToGlobal(tgContext, rPhiOnDisc, ignoredMomentum,
                                   returnedPosition);
  CHECK_CLOSE_ABS(returnedPosition, expectedPosition, 1e-6);
  //
  discSurfaceObject->localToGlobal(tgContext, rPhiNotInSector, ignoredMomentum,
                                   returnedPosition);
  Vector3D expectedNonPosition{-1.2, 0, 0};
  CHECK_CLOSE_ABS(returnedPosition, expectedNonPosition, 1e-6);
  //
  /// Test globalToLocal
  Vector2D returnedLocalPosition{33., 44.};
  Vector2D expectedLocalPosition{1.2, 0.0};
  BOOST_CHECK(discSurfaceObject->globalToLocal(tgContext, point3DOnSurface,
                                               ignoredMomentum,
                                               returnedLocalPosition));  // pass
  CHECK_CLOSE_ABS(returnedLocalPosition, expectedLocalPosition, 1e-6);

  // Global to local does not check inside bounds
  BOOST_CHECK(discSurfaceObject->globalToLocal(
      tgContext, point3DNotInSector, ignoredMomentum, returnedLocalPosition));
  //
  Vector3D pointOutsideR{0.0, 100., 0};
  BOOST_CHECK(discSurfaceObject->globalToLocal(
      tgContext, pointOutsideR, ignoredMomentum, returnedLocalPosition));
  //
  /// Test localPolarToCartesian
  Vector2D rPhi1_1{std::sqrt(2.), M_PI / 4.};
  Vector2D cartesian1_1{1., 1.};
  CHECK_CLOSE_REL(discSurfaceObject->localPolarToCartesian(rPhi1_1),
                  cartesian1_1, 1e-6);
  //
  /// Test localCartesianToPolar
  CHECK_CLOSE_REL(discSurfaceObject->localCartesianToPolar(cartesian1_1),
                  rPhi1_1, 1e-6);
  //
  /// Test localPolarToLocalCartesian
  CHECK_CLOSE_REL(discSurfaceObject->localPolarToLocalCartesian(rPhi1_1),
                  cartesian1_1, 1e-6);
  //
  /// Test localCartesianToGlobal
  Vector3D cartesian3D1_1{1., 1., 0.};
  CHECK_CLOSE_ABS(
      discSurfaceObject->localCartesianToGlobal(tgContext, cartesian1_1),
      cartesian3D1_1, 1e-6);
  //
  /// Test globalToLocalCartesian
  CHECK_CLOSE_REL(
      discSurfaceObject->globalToLocalCartesian(tgContext, cartesian3D1_1),
      cartesian1_1, 1e-6);
  //
  /// Test pathCorrection
  double projected3DMomentum = std::sqrt(3.) * 1.e6;
  Vector3D momentum{projected3DMomentum, projected3DMomentum,
                    projected3DMomentum};
  Vector3D ignoredPosition{1.1, 2.2, 3.3};
  CHECK_CLOSE_REL(discSurfaceObject->pathCorrection(tgContext, ignoredPosition,
                                                    momentum.normalized()),
                  std::sqrt(3), 0.01);
  //
  /// intersectionEstimate
  Vector3D globalPosition{1.2, 0.0, -10.};
  Vector3D direction{0., 0., 1.};  // must be normalised
  Vector3D expected{1.2, 0.0, 0.0};
  // intersect is a struct of (Vector3D) position, pathLength, distance and
  // (bool) valid
  auto intersect = discSurfaceObject->intersectionEstimate(
      tgContext, globalPosition, direction, false);
  Intersection expectedIntersect{Vector3D{1.2, 0., 0.}, 10.,
                                 Intersection::Status::reachable};
  BOOST_CHECK(bool(intersect));
  CHECK_CLOSE_ABS(intersect.position, expectedIntersect.position, 1e-9);
  CHECK_CLOSE_ABS(intersect.pathLength, expectedIntersect.pathLength, 1e-9);
  //
  /// Test name
  boost::test_tools::output_test_stream nameOuput;
  nameOuput << discSurfaceObject->name();
  BOOST_CHECK(nameOuput.is_equal("Acts::DiscSurface"));
}
//
/// Unit test for testing DiscSurface assignment and equality
BOOST_AUTO_TEST_CASE(DiscSurfaceAssignment) {
  Vector3D origin3D{0, 0, 0};
  std::shared_ptr<const Transform3D> pTransform;  // nullptr
  double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
  auto discSurfaceObject =
      Surface::makeShared<DiscSurface>(pTransform, rMin, rMax, halfPhiSector);
  auto assignedDisc = Surface::makeShared<DiscSurface>(nullptr, 2.2, 4.4, 0.07);
  //
  BOOST_CHECK_NO_THROW(*assignedDisc = *discSurfaceObject);
  BOOST_CHECK((*assignedDisc) == (*discSurfaceObject));
}

/// Unit test for testing DiscSurface assignment and equality
BOOST_AUTO_TEST_CASE(DiscSurfaceExtent) {
  double rMin(1.0), rMax(5.0);

  auto pDisc = Surface::makeShared<DiscSurface>(nullptr, 0., rMax);
  auto pDiscExtent = pDisc->polyhedronRepresentation(tgContext, 1).extent();

  CHECK_CLOSE_ABS(0., pDiscExtent.min(binZ), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(0., pDiscExtent.max(binZ), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(0., pDiscExtent.min(binR), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(rMax, pDiscExtent.max(binR), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(-rMax, pDiscExtent.min(binX), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(rMax, pDiscExtent.max(binX), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(-rMax, pDiscExtent.min(binY), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(rMax, pDiscExtent.max(binY), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(-M_PI, pDiscExtent.min(binPhi), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(M_PI, pDiscExtent.max(binPhi), s_onSurfaceTolerance);

  auto pRing = Surface::makeShared<DiscSurface>(nullptr, rMin, rMax);
  auto pRingExtent = pRing->polyhedronRepresentation(tgContext, 1).extent();

  CHECK_CLOSE_ABS(0., pRingExtent.min(binZ), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(0., pRingExtent.max(binZ), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(rMin, pRingExtent.min(binR), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(rMax, pRingExtent.max(binR), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(-rMax, pRingExtent.min(binX), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(rMax, pRingExtent.max(binX), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(-rMax, pRingExtent.min(binY), s_onSurfaceTolerance);
  CHECK_CLOSE_ABS(rMax, pRingExtent.max(binY), s_onSurfaceTolerance);
}

/// Unit test for testing DiscSurface alignment derivatives
BOOST_AUTO_TEST_CASE(DiscSurfaceAlignment) {
  Translation3D translation{0., 1., 2.};
  Transform3D transform(translation);
  auto pTransform = std::make_shared<const Transform3D>(translation);
  double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
  auto discSurfaceObject =
      Surface::makeShared<DiscSurface>(pTransform, rMin, rMax, halfPhiSector);

  const auto& rotation = pTransform->rotation();
  // The local frame z axis
  const Vector3D localZAxis = rotation.col(2);
  // Check the local z axis is aligned to global z axis
  CHECK_CLOSE_ABS(localZAxis, Vector3D(0., 0., 1.), 1e-15);

  /// Define the track (global) position and direction
  Vector3D globalPosition{0, 4, 2};
  Vector3D momentum{0, 0, 1};
  Vector3D direction = momentum.normalized();

  // Call the function to calculate the derivative of local frame axes w.r.t its
  // rotation
  const auto& [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);

  // (a) Test the derivative of path length w.r.t. alignment parameters
  const AlignmentRowVector& alignToPath =
      discSurfaceObject->alignmentToPathDerivative(tgContext, rotToLocalZAxis,
                                                   globalPosition, direction);
  // The expected results
  AlignmentRowVector expAlignToPath = AlignmentRowVector::Zero();
  expAlignToPath << 0, 0, 1, 3, 0, 0;
  // Check if the calculated derivative is as expected
  CHECK_CLOSE_ABS(alignToPath, expAlignToPath, 1e-10);

  // (b) Test the derivative of bound track parameters local position w.r.t.
  // position in local 3D Cartesian coordinates
  const auto& loc3DToLocBound =
      discSurfaceObject->localCartesianToBoundLocalDerivative(tgContext,
                                                              globalPosition);
  // Check if the result is as expected
  LocalCartesianToBoundLocalMatrix expLoc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Zero();
  expLoc3DToLocBound << 0, 1, 0, -1.0 / 3, 0, 0;
  CHECK_CLOSE_ABS(loc3DToLocBound, expLoc3DToLocBound, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
