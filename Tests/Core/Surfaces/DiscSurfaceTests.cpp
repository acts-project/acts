// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Disc Surface Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidalBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace utf = boost::unit_test;
namespace tt  = boost::test_tools;

namespace Acts {

namespace Test {
  // using boost::test_tools::output_test_stream;

  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit tests for creating DiscSurface object
  BOOST_AUTO_TEST_CASE(DiscSurface_constructors_test)
  {
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
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    BOOST_CHECK_NO_THROW(Surface::makeShared<DiscSurface>(
        pTransform, rMin, rMax, halfPhiSector));
    //
    /// Copy constructed DiscSurface
    auto anotherDiscSurface = Surface::makeShared<DiscSurface>(
        pTransform, rMin, rMax, halfPhiSector);
    // N.B. Just using
    // BOOST_CHECK_NO_THROW(Surface::makeShared<DiscSurface>(anotherDiscSurface))
    // tries to call
    // the (deleted) default constructor.
    auto copiedSurface = Surface::makeShared<DiscSurface>(*anotherDiscSurface);
    BOOST_TEST_MESSAGE("Copy constructed DiscSurface ok");
    //
    /// Copied and transformed DiscSurface
    BOOST_CHECK_NO_THROW(
        Surface::makeShared<DiscSurface>(*anotherDiscSurface, *pTransform));

    /// Construct with nullptr bounds
    DetectorElementStub detElem;
    BOOST_CHECK_THROW(auto nullBounds
                      = Surface::makeShared<DiscSurface>(nullptr, detElem),
                      AssertionFailureException);
  }

  /// Unit tests of all named methods
  BOOST_AUTO_TEST_CASE(DiscSurface_properties_test, *utf::expected_failures(2))
  {
    Vector3D                           origin3D{0, 0, 0};
    std::shared_ptr<const Transform3D> pTransform;  // nullptr
    double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
    auto   discSurfaceObject = Surface::makeShared<DiscSurface>(
        pTransform, rMin, rMax, halfPhiSector);
    //
    /// Test type
    BOOST_CHECK_EQUAL(discSurfaceObject->type(), Surface::Disc);
    //
    /// Test normal, no local position specified
    Vector3D zAxis{0, 0, 1};
    BOOST_CHECK_EQUAL(discSurfaceObject->normal(), zAxis);
    //
    /// Test normal, local position specified
    Vector2D lpos(2.0, 0.05);
    BOOST_CHECK_EQUAL(discSurfaceObject->normal(lpos), zAxis);
    //
    /// Test binningPosition
    // auto binningPosition=
    // discSurfaceObject.binningPosition(BinningValue::binRPhi );
    // std::cout<<binningPosition<<std::endl;
    BOOST_CHECK_EQUAL(discSurfaceObject->binningPosition(BinningValue::binRPhi),
                      origin3D);
    //
    /// Test bounds
    BOOST_CHECK_EQUAL(discSurfaceObject->bounds().type(), SurfaceBounds::Disc);
    //
    Vector3D ignoredMomentum{0., 0., 0.};
    /// Test isOnSurface()
    Vector3D point3DNotInSector{0.0, 1.2, 0};
    Vector3D point3DOnSurface{1.2, 0.0, 0};
    BOOST_CHECK(!discSurfaceObject->isOnSurface(
        point3DNotInSector, ignoredMomentum, true));  // passes
    BOOST_CHECK(discSurfaceObject->isOnSurface(
        point3DOnSurface, ignoredMomentum, true));  // passes
    //
    /// Test localToGlobal
    Vector3D returnedPosition{10.9, 8.7, 6.5};
    Vector3D expectedPosition{1.2, 0, 0};
    Vector2D rPhiOnDisc{1.2, 0.0};
    Vector2D rPhiNotInSector{1.2, M_PI};  // outside sector at Phi=0, +/- pi/8
    discSurfaceObject->localToGlobal(
        rPhiOnDisc, ignoredMomentum, returnedPosition);
    CHECK_CLOSE_ABS(returnedPosition, expectedPosition, 1e-6);
    //
    discSurfaceObject->localToGlobal(
        rPhiNotInSector, ignoredMomentum, returnedPosition);
    Vector3D expectedNonPosition{-1.2, 0, 0};
    CHECK_CLOSE_ABS(returnedPosition, expectedNonPosition, 1e-6);
    //
    /// Test globalToLocal
    Vector2D returnedLocalPosition{33., 44.};
    Vector2D expectedLocalPosition{1.2, 0.0};
    BOOST_CHECK(discSurfaceObject->globalToLocal(
        point3DOnSurface, ignoredMomentum, returnedLocalPosition));  // pass
    CHECK_CLOSE_ABS(returnedLocalPosition, expectedLocalPosition, 1e-6);
    //
    BOOST_CHECK(!discSurfaceObject->globalToLocal(
        point3DNotInSector,
        ignoredMomentum,
        returnedLocalPosition));  // test fails
    //
    Vector3D pointOutsideRadius{0.0, 100., 0};
    BOOST_CHECK(
        !discSurfaceObject->globalToLocal(pointOutsideRadius,
                                          ignoredMomentum,
                                          returnedLocalPosition));  // fails
    //
    /// Test localPolarToCartesian
    Vector2D rPhi1_1{std::sqrt(2.), M_PI / 4.};
    Vector2D cartesian1_1{1., 1.};
    CHECK_CLOSE_REL(
        discSurfaceObject->localPolarToCartesian(rPhi1_1), cartesian1_1, 1e-6);
    //
    /// Test localCartesianToPolar
    CHECK_CLOSE_REL(
        discSurfaceObject->localCartesianToPolar(cartesian1_1), rPhi1_1, 1e-6);
    //
    /// Test localPolarToLocalCartesian
    CHECK_CLOSE_REL(discSurfaceObject->localPolarToLocalCartesian(rPhi1_1),
                    cartesian1_1,
                    1e-6);
    //
    /// Test localCartesianToGlobal
    Vector3D cartesian3D1_1{1., 1., 0.};
    CHECK_CLOSE_ABS(discSurfaceObject->localCartesianToGlobal(cartesian1_1),
                    cartesian3D1_1,
                    1e-6);
    //
    /// Test globalToLocalCartesian
    CHECK_CLOSE_REL(discSurfaceObject->globalToLocalCartesian(cartesian3D1_1),
                    cartesian1_1,
                    1e-6);
    //
    /// Test pathCorrection
    double   projected3DMomentum = std::sqrt(3.) * 1.e6;
    Vector3D momentum{
        projected3DMomentum, projected3DMomentum, projected3DMomentum};
    Vector3D ignoredPosition{1.1, 2.2, 3.3};
    CHECK_CLOSE_REL(
        discSurfaceObject->pathCorrection(ignoredPosition, momentum),
        std::sqrt(3),
        0.01);
    //
    /// intersectionEstimate
    Vector3D globalPosition{1.2, 0.0, -10.};
    Vector3D direction{0., 0., 1.};  // must be normalised
    Vector3D expected{1.2, 0.0, 0.0};
    // intersect is a struct of (Vector3D) position, pathLength, distance and
    // (bool) valid
    auto intersect
        = discSurfaceObject->intersectionEstimate(globalPosition, direction);
    Intersection expectedIntersect{Vector3D{1.2, 0., 0.}, 10., true, 0.0};
    BOOST_CHECK(intersect.valid);
    CHECK_CLOSE_ABS(intersect.position, expectedIntersect.position, 1e-9);
    CHECK_CLOSE_ABS(intersect.pathLength, expectedIntersect.pathLength, 1e-9);
    CHECK_CLOSE_ABS(intersect.distance, expectedIntersect.distance, 1e-9);
    //
    /// Test name
    boost::test_tools::output_test_stream nameOuput;
    nameOuput << discSurfaceObject->name();
    BOOST_CHECK(nameOuput.is_equal("Acts::DiscSurface"));
  }
  //
  /// Unit test for testing DiscSurface assignment and equality
  BOOST_AUTO_TEST_CASE(DiscSurface_assignment_test)
  {
    Vector3D                           origin3D{0, 0, 0};
    std::shared_ptr<const Transform3D> pTransform;  // nullptr
    double rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
    auto   discSurfaceObject = Surface::makeShared<DiscSurface>(
        pTransform, rMin, rMax, halfPhiSector);
    auto assignedDisc
        = Surface::makeShared<DiscSurface>(nullptr, 2.2, 4.4, 0.07);
    //
    BOOST_CHECK_NO_THROW(*assignedDisc = *discSurfaceObject);
    BOOST_CHECK_EQUAL(*assignedDisc, *discSurfaceObject);
  }

  BOOST_AUTO_TEST_CASE(DiscSurface_toVariantData)
  {
    double minR = 1, maxR = 4, avgPhi = M_PI / 3, phiSec = M_PI;
    auto   rbounds
        = std::make_shared<const RadialBounds>(minR, maxR, avgPhi, phiSec);

    Transform3D rot(AngleAxis3D(M_PI / 4., Vector3D::UnitZ()));
    auto        trf
        = std::make_shared<const Transform3D>(Translation3D(0, 0, 2) * rot);

    auto rdisc = Surface::makeShared<DiscSurface>(trf, rbounds);

    variant_data var_rdisc = rdisc->toVariantData();
    std::cout << var_rdisc << std::endl;

    auto rdisc2 = Surface::makeShared<DiscSurface>(var_rdisc);
    CHECK_CLOSE_OR_SMALL(rdisc->transform(), rdisc2->transform(), 1e-6, 1e-9);
    CHECK_CLOSE_OR_SMALL(rdisc2->transform(), *trf, 1e-6, 1e-9);

    const RadialBounds* rbounds_act
        = dynamic_cast<const RadialBounds*>(&rdisc2->bounds());
    BOOST_CHECK_EQUAL(rbounds->rMin(), rbounds_act->rMin());
    BOOST_CHECK_EQUAL(rbounds->rMax(), rbounds_act->rMax());
    BOOST_CHECK_EQUAL(rbounds->averagePhi(), rbounds_act->averagePhi());
    BOOST_CHECK_EQUAL(rbounds->halfPhiSector(), rbounds_act->halfPhiSector());

    double rMin = 1, rMax = 5, minHalfX = 2, maxHalfX = 4, stereo = M_PI / 8.;
    auto   dtbounds = std::make_shared<const DiscTrapezoidalBounds>(
        minHalfX, maxHalfX, rMin, rMax, avgPhi, stereo);

    auto         dtdisc     = Surface::makeShared<DiscSurface>(trf, dtbounds);
    variant_data var_dtdisc = dtdisc->toVariantData();
    std::cout << var_dtdisc;

    auto dtdisc2 = Surface::makeShared<DiscSurface>(var_dtdisc);

    CHECK_CLOSE_OR_SMALL(dtdisc->transform(), dtdisc2->transform(), 1e-6, 1e-9);
    CHECK_CLOSE_OR_SMALL(dtdisc2->transform(), *trf, 1e-6, 1e-9);

    const DiscTrapezoidalBounds* dtbounds_act
        = dynamic_cast<const DiscTrapezoidalBounds*>(&dtdisc2->bounds());
    BOOST_CHECK_EQUAL(dtbounds->rMin(), dtbounds_act->rMin());
    BOOST_CHECK_EQUAL(dtbounds->rMax(), dtbounds_act->rMax());
    BOOST_CHECK_EQUAL(dtbounds->minHalflengthX(),
                      dtbounds_act->minHalflengthX());
    BOOST_CHECK_EQUAL(dtbounds->maxHalflengthX(),
                      dtbounds_act->maxHalflengthX());
    BOOST_CHECK_EQUAL(dtbounds->stereo(), dtbounds_act->stereo());
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
