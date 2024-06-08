// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
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
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit tests for creating DiscSurface object
BOOST_AUTO_TEST_CASE(DiscSurfaceConstruction) {
  // default constructor is deleted
  // scaffolding...
  double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
  //
  /// Test DiscSurface constructor with default halfPhiSector
  BOOST_CHECK_NO_THROW(
      Surface::makeShared<DiscSurface>(Transform3::Identity(), rMin, rMax));
  //
  /// Test DiscSurface constructor with a transform specified
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
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
      tgContext, *anotherDiscSurface, pTransform));

  /// Construct with nullptr bounds
  DetectorElementStub detElem;
  BOOST_CHECK_THROW(
      auto nullBounds = Surface::makeShared<DiscSurface>(nullptr, detElem),
      AssertionFailureException);
}

/// Unit tests of all named methods
BOOST_AUTO_TEST_CASE(DiscSurfaceProperties) {
  Vector3 origin3D{0, 0, 0};
  double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
  auto discSurfaceObject = Surface::makeShared<DiscSurface>(
      Transform3::Identity(), rMin, rMax, halfPhiSector);
  //
  /// Test type
  BOOST_CHECK_EQUAL(discSurfaceObject->type(), Surface::Disc);
  //
  /// Test normal, no local position specified
  Vector3 zAxis{0, 0, 1};
  BOOST_CHECK_EQUAL(discSurfaceObject->normal(tgContext), zAxis);
  //
  /// Test normal, local position specified
  Vector2 lpos(2.0, 0.05);
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
  Vector3 ignoredMomentum{0., 0., 0.};
  /// Test isOnSurface()
  Vector3 point3DNotInSector{0.0, 1.2, 0};
  Vector3 point3DOnSurface{1.2, 0.0, 0};
  BOOST_CHECK(!discSurfaceObject->isOnSurface(tgContext, point3DNotInSector,
                                              ignoredMomentum,
                                              BoundaryCheck(true)));  // passes
  BOOST_CHECK(discSurfaceObject->isOnSurface(tgContext, point3DOnSurface,
                                             ignoredMomentum,
                                             BoundaryCheck(true)));  // passes
  //
  /// Test localToGlobal
  Vector3 returnedPosition{10.9, 8.7, 6.5};
  Vector3 expectedPosition{1.2, 0, 0};
  Vector2 rPhiOnDisc{1.2, 0.0};
  Vector2 rPhiNotInSector{1.2, M_PI};  // outside sector at Phi=0, +/- pi/8
  returnedPosition =
      discSurfaceObject->localToGlobal(tgContext, rPhiOnDisc, ignoredMomentum);
  CHECK_CLOSE_ABS(returnedPosition, expectedPosition, 1e-6);
  //
  returnedPosition = discSurfaceObject->localToGlobal(
      tgContext, rPhiNotInSector, ignoredMomentum);
  Vector3 expectedNonPosition{-1.2, 0, 0};
  CHECK_CLOSE_ABS(returnedPosition, expectedNonPosition, 1e-6);
  //
  /// Test globalToLocal
  Vector2 returnedLocalPosition{33., 44.};
  Vector2 expectedLocalPosition{1.2, 0.0};
  returnedLocalPosition =
      discSurfaceObject
          ->globalToLocal(tgContext, point3DOnSurface, ignoredMomentum)
          .value();
  CHECK_CLOSE_ABS(returnedLocalPosition, expectedLocalPosition, 1e-6);

  // Global to local does not check inside bounds
  returnedLocalPosition =
      discSurfaceObject
          ->globalToLocal(tgContext, point3DNotInSector, ignoredMomentum)
          .value();
  //
  Vector3 pointOutsideR{0.0, 100., 0};
  returnedLocalPosition =
      discSurfaceObject
          ->globalToLocal(tgContext, pointOutsideR, ignoredMomentum)
          .value();
  //
  /// Test localPolarToCartesian
  Vector2 rPhi1_1{std::sqrt(2.), M_PI / 4.};
  Vector2 cartesian1_1{1., 1.};
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
  Vector3 cartesian3D1_1{1., 1., 0.};
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
  Vector3 momentum{projected3DMomentum, projected3DMomentum,
                   projected3DMomentum};
  Vector3 ignoredPosition = discSurfaceObject->center(tgContext);
  CHECK_CLOSE_REL(discSurfaceObject->pathCorrection(tgContext, ignoredPosition,
                                                    momentum.normalized()),
                  std::sqrt(3), 0.01);
  //
  /// intersection test
  Vector3 globalPosition{1.2, 0.0, -10.};
  Vector3 direction{0., 0., 1.};  // must be normalised
  Vector3 expected{1.2, 0.0, 0.0};

  // intersect is a struct of (Vector3) position, pathLength, distance and
  // (bool) valid, it's contained in a Surface intersection
  auto sfIntersection = discSurfaceObject
                            ->intersect(tgContext, globalPosition, direction,
                                        BoundaryCheck(false))
                            .closest();
  Intersection3D expectedIntersect{Vector3{1.2, 0., 0.}, 10.,
                                   Intersection3D::Status::reachable};
  BOOST_CHECK(bool(sfIntersection));
  CHECK_CLOSE_ABS(sfIntersection.position(), expectedIntersect.position(),
                  1e-9);
  CHECK_CLOSE_ABS(sfIntersection.pathLength(), expectedIntersect.pathLength(),
                  1e-9);
  BOOST_CHECK_EQUAL(sfIntersection.object(), discSurfaceObject.get());

  //
  /// Test name
  boost::test_tools::output_test_stream nameOuput;
  nameOuput << discSurfaceObject->name();
  BOOST_CHECK(nameOuput.is_equal("Acts::DiscSurface"));
}
//
/// Unit test for testing DiscSurface assignment and equality
BOOST_AUTO_TEST_CASE(DiscSurfaceAssignment) {
  Vector3 origin3D{0, 0, 0};
  double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
  auto discSurfaceObject = Surface::makeShared<DiscSurface>(
      Transform3::Identity(), rMin, rMax, halfPhiSector);
  auto assignedDisc =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 2.2, 4.4, 0.07);
  //
  BOOST_CHECK_NO_THROW(*assignedDisc = *discSurfaceObject);
  BOOST_CHECK((*assignedDisc) == (*discSurfaceObject));
}

/// Unit test for testing DiscSurface assignment and equality
BOOST_AUTO_TEST_CASE(DiscSurfaceExtent) {
  double rMin(1.0), rMax(5.0);

  auto pDisc =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 0., rMax);
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

  auto pRing =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), rMin, rMax);
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
  Translation3 translation{0., 1., 2.};
  Transform3 transform(translation);
  double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
  auto discSurfaceObject =
      Surface::makeShared<DiscSurface>(transform, rMin, rMax, halfPhiSector);

  const auto& rotation = transform.rotation();
  // The local frame z axis
  const Vector3 localZAxis = rotation.col(2);
  // Check the local z axis is aligned to global z axis
  CHECK_CLOSE_ABS(localZAxis, Vector3(0., 0., 1.), 1e-15);

  // Define the track (global) position and direction
  Vector3 globalPosition{0, 4, 2};
  Vector3 momentum{0, 0, 1};
  Vector3 direction = momentum.normalized();

  // (a) Test the derivative of path length w.r.t. alignment parameters
  const AlignmentToPathMatrix& alignToPath =
      discSurfaceObject->alignmentToPathDerivative(tgContext, globalPosition,
                                                   direction);
  // The expected results
  AlignmentToPathMatrix expAlignToPath = AlignmentToPathMatrix::Zero();
  expAlignToPath << 0, 0, 1, 3, 0, 0;
  // Check if the calculated derivative is as expected
  CHECK_CLOSE_ABS(alignToPath, expAlignToPath, 1e-10);

  // (b) Test the derivative of bound track parameters local position w.r.t.
  // position in local 3D Cartesian coordinates
  const auto& loc3DToLocBound =
      discSurfaceObject->localCartesianToBoundLocalDerivative(tgContext,
                                                              globalPosition);
  // Check if the result is as expected
  ActsMatrix<2, 3> expLoc3DToLocBound = ActsMatrix<2, 3>::Zero();
  expLoc3DToLocBound << 0, 1, 0, -1.0 / 3, 0, 0;
  CHECK_CLOSE_ABS(loc3DToLocBound, expLoc3DToLocBound, 1e-10);
}

BOOST_AUTO_TEST_CASE(DiscSurfaceBinningPosition) {
  using namespace Acts::UnitLiterals;
  Vector3 s{5_mm, 7_mm, 10_cm};
  Transform3 trf;
  trf = Translation3(s) * AngleAxis3{0.5, Vector3::UnitZ()};

  double minR = 300;
  double maxR = 330;

  {
    // Radial Bounds
    auto bounds = std::make_shared<RadialBounds>(minR, maxR, M_PI / 8, 0.1);
    auto disc = Acts::Surface::makeShared<Acts::DiscSurface>(trf, bounds);

    Vector3 bp = disc->binningPosition(tgContext, binR);
    double r = (bounds->rMax() + bounds->rMin()) / 2.0;
    double phi = bounds->get(RadialBounds::eAveragePhi);
    Vector3 exp = Vector3{r * std::cos(phi), r * std::sin(phi), 0};
    exp = trf * exp;

    BOOST_CHECK_EQUAL(bp, exp);
    BOOST_CHECK_EQUAL(disc->binningPositionValue(tgContext, binR),
                      VectorHelpers::perp(exp));

    bp = disc->binningPosition(tgContext, binPhi);
    BOOST_CHECK_EQUAL(bp, exp);
    BOOST_CHECK_EQUAL(disc->binningPositionValue(tgContext, binPhi),
                      VectorHelpers::phi(exp));

    for (auto b : {binX, binY, binZ, binEta, binRPhi, binH, binMag}) {
      BOOST_TEST_CONTEXT("binValue: " << b) {
        BOOST_CHECK_EQUAL(disc->binningPosition(tgContext, b),
                          disc->center(tgContext));
      }
    }
  }

  {
    // Annulus Bounds
    double minPhiRel = -0.3;
    double maxPhiRel = 0.2;
    Vector2 origin{5_mm, 5_mm};
    auto bounds = std::make_shared<AnnulusBounds>(minR, maxR, minPhiRel,
                                                  maxPhiRel, origin);

    auto disc = Acts::Surface::makeShared<Acts::DiscSurface>(trf, bounds);

    Vector3 bp = disc->binningPosition(tgContext, binR);
    double r = (bounds->rMax() + bounds->rMin()) / 2.0;
    double phi = bounds->get(AnnulusBounds::eAveragePhi);
    Vector3 exp = Vector3{r * std::cos(phi), r * std::sin(phi), 0};
    exp = trf * exp;

    BOOST_CHECK_EQUAL(bp, exp);

    bp = disc->binningPosition(tgContext, binPhi);
    BOOST_CHECK_EQUAL(bp, exp);

    for (auto b : {binX, binY, binZ, binEta, binRPhi, binH, binMag}) {
      BOOST_TEST_CONTEXT("binValue: " << b) {
        BOOST_CHECK_EQUAL(disc->binningPosition(tgContext, b),
                          disc->center(tgContext));
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
