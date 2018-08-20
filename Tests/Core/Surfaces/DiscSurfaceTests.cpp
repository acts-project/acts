// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Disc Surface Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/DiscTrapezoidalBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"
//
#include "DetectorElementStub.hpp"
//
#include <limits>

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
    BOOST_CHECK_NO_THROW(DiscSurface(nullptr, rMin, rMax, halfPhiSector));
    //
    /// Test DiscSurface constructor with default halfPhiSector
    BOOST_CHECK_NO_THROW(DiscSurface(nullptr, rMin, rMax));
    //
    /// Test DiscSurface constructor with a transform specified
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    BOOST_CHECK_NO_THROW(DiscSurface(pTransform, rMin, rMax, halfPhiSector));
    //
    /// Copy constructed DiscSurface
    DiscSurface anotherDiscSurface(pTransform, rMin, rMax, halfPhiSector);
    // N.B. Just using BOOST_CHECK_NO_THROW(DiscSurface(anotherDiscSurface))
    // tries to call
    // the (deleted) default constructor.
    DiscSurface copiedDiscSurface(anotherDiscSurface);
    BOOST_TEST_MESSAGE("Copy constructed DiscSurface ok");
    //
    /// Copied and transformed DiscSurface
    BOOST_CHECK_NO_THROW(DiscSurface(anotherDiscSurface, *pTransform));

    /// Construct with nullptr bounds
    DetectorElementStub detElem;
    BOOST_CHECK_THROW(DiscSurface nullBounds(nullptr, detElem),
                      AssertionFailureException);
  }

  /// Unit tests of all named methods
  BOOST_AUTO_TEST_CASE(DiscSurface_properties_test, *utf::expected_failures(2))
  {
    Vector3D                           origin3D{0, 0, 0};
    std::shared_ptr<const Transform3D> pTransform;  // nullptr
    double      rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
    DiscSurface discSurfaceObject(pTransform, rMin, rMax, halfPhiSector);
    //
    /// Test type
    BOOST_CHECK(discSurfaceObject.type() == Surface::Disc);
    //
    /// Test normal, no local position specified
    Vector3D zAxis{0, 0, 1};
    BOOST_CHECK(discSurfaceObject.normal() == zAxis);
    //
    /// Test normal, local position specified
    Vector2D lpos(2.0, 0.05);
    BOOST_CHECK(discSurfaceObject.normal(lpos) == zAxis);
    //
    /// Test binningPosition
    // auto binningPosition=
    // discSurfaceObject.binningPosition(BinningValue::binRPhi );
    // std::cout<<binningPosition<<std::endl;
    BOOST_CHECK(discSurfaceObject.binningPosition(BinningValue::binRPhi)
                == origin3D);
    //
    /// Test bounds
    BOOST_TEST(discSurfaceObject.bounds().type() = SurfaceBounds::Disc);
    //
    /// Test isOnSurface()
    Vector3D point3DNotInSector{0.0, 1.2, 0};
    Vector3D point3DOnSurface{1.2, 0.0, 0};
    BOOST_TEST(discSurfaceObject.isOnSurface(point3DNotInSector, true)
               == false);  // passes
    BOOST_TEST(discSurfaceObject.isOnSurface(point3DOnSurface, true)
               == true);  // passes
    //
    /// Test localToGlobal
    Vector3D returnedPosition{10.9, 8.7, 6.5};
    Vector3D expectedPosition{1.2, 0, 0};
    Vector2D rPhiOnDisc{1.2, 0.0};
    Vector2D rPhiNotInSector{1.2, M_PI};  // outside sector at Phi=0, +/- pi/8
    Vector3D ignoredMomentum{0., 0., 0.};
    discSurfaceObject.localToGlobal(
        rPhiOnDisc, ignoredMomentum, returnedPosition);
    BOOST_TEST(returnedPosition.isApprox(expectedPosition),
               "LocalToGlobal for rPhiOnDisc");
    //
    discSurfaceObject.localToGlobal(
        rPhiNotInSector, ignoredMomentum, returnedPosition);
    Vector3D expectedNonPosition{-1.2, 0, 0};
    BOOST_TEST(returnedPosition.isApprox(expectedNonPosition));
    //
    /// Test globalToLocal
    Vector2D returnedLocalPosition{33., 44.};
    Vector2D expectedLocalPosition{1.2, 0.0};
    BOOST_TEST(discSurfaceObject.globalToLocal(
        point3DOnSurface, ignoredMomentum, returnedLocalPosition));  // pass
    BOOST_TEST(returnedLocalPosition.isApprox(expectedLocalPosition));
    //
    BOOST_TEST(discSurfaceObject.globalToLocal(
                   point3DNotInSector, ignoredMomentum, returnedLocalPosition)
               == false);  // test fails
    //
    Vector3D pointOutsideRadius{0.0, 100., 0};
    BOOST_TEST(discSurfaceObject.globalToLocal(
                   pointOutsideRadius, ignoredMomentum, returnedLocalPosition)
               == false);  // fails
    //
    /// Test localPolarToCartesian
    Vector2D rPhi1_1{std::sqrt(2.), M_PI / 4.};
    Vector2D cartesian1_1{1., 1.};
    BOOST_TEST(discSurfaceObject.localPolarToCartesian(rPhi1_1).isApprox(
        cartesian1_1));
    //
    /// Test localCartesianToPolar
    BOOST_TEST(discSurfaceObject.localCartesianToPolar(cartesian1_1)
                   .isApprox(rPhi1_1));
    //
    /// Test localPolarToLocalCartesian
    BOOST_TEST(discSurfaceObject.localPolarToLocalCartesian(rPhi1_1).isApprox(
        cartesian1_1));
    //
    /// Test localCartesianToGlobal
    Vector3D cartesian3D1_1{1., 1., 0.};
    BOOST_TEST(discSurfaceObject.localCartesianToGlobal(cartesian1_1)
                   .isApprox(cartesian3D1_1));
    //
    /// Test globalToLocalCartesian
    BOOST_TEST(discSurfaceObject.globalToLocalCartesian(cartesian3D1_1)
                   .isApprox(cartesian1_1));
    //
    /// Test pathCorrection
    double   projected3DMomentum = std::sqrt(3.) * 1.e6;
    Vector3D momentum{
        projected3DMomentum, projected3DMomentum, projected3DMomentum};
    Vector3D ignoredPosition{1.1, 2.2, 3.3};
    BOOST_TEST(discSurfaceObject.pathCorrection(ignoredPosition, momentum)
                   == std::sqrt(3),
               tt::tolerance(0.01));
    //
    /// intersectionEstimate
    Vector3D globalPosition{1.2, 0.0, -10.};
    Vector3D direction{0., 0., 1.};  // must be normalised
    Vector3D expected{1.2, 0.0, 0.0};
    // intersect is a struct of (Vector3D) position, pathLength, distance and
    // (bool) valid
    auto intersect
        = discSurfaceObject.intersectionEstimate(globalPosition, direction);
    Intersection expectedIntersect{Vector3D{1.2, 0., 0.}, 10., true, 0.0};
    BOOST_TEST(intersect.valid);
    BOOST_TEST(intersect.position.isApprox(expectedIntersect.position));
    BOOST_TEST(intersect.pathLength == expectedIntersect.pathLength);
    BOOST_TEST(intersect.distance == expectedIntersect.distance);
    //
    /// Test name
    boost::test_tools::output_test_stream nameOuput;
    nameOuput << discSurfaceObject.name();
    BOOST_TEST(nameOuput.is_equal("Acts::DiscSurface"));
  }
  //
  /// Unit test for testing DiscSurface assignment and equality
  BOOST_AUTO_TEST_CASE(DiscSurface_assignment_test)
  {
    Vector3D                           origin3D{0, 0, 0};
    std::shared_ptr<const Transform3D> pTransform;  // nullptr
    double      rMin(1.0), rMax(5.0), halfPhiSector(M_PI / 8.);
    DiscSurface discSurfaceObject(pTransform, rMin, rMax, halfPhiSector);
    DiscSurface assignedDisc(nullptr, 2.2, 4.4, 0.07);
    //
    BOOST_CHECK_NO_THROW(assignedDisc = discSurfaceObject);
    BOOST_CHECK(assignedDisc == discSurfaceObject);
  }

  BOOST_AUTO_TEST_CASE(DiscSurface_toVariantData)
  {
    double minR = 1, maxR = 4, avgPhi = M_PI / 3, phiSec = M_PI;
    auto   rbounds
        = std::make_shared<const RadialBounds>(minR, maxR, avgPhi, phiSec);

    Transform3D rot(AngleAxis3D(M_PI / 4., Vector3D::UnitZ()));
    auto        trf
        = std::make_shared<const Transform3D>(Translation3D(0, 0, 2) * rot);

    DiscSurface rdisc(trf, rbounds);

    variant_data var_rdisc = rdisc.toVariantData();
    std::cout << var_rdisc << std::endl;

    DiscSurface rdisc2(var_rdisc);
    BOOST_TEST(rdisc.transform().isApprox(rdisc2.transform()));
    BOOST_TEST(rdisc2.transform().isApprox(*trf));

    const RadialBounds* rbounds_act
        = dynamic_cast<const RadialBounds*>(&rdisc2.bounds());
    BOOST_TEST(rbounds->rMin() == rbounds_act->rMin());
    BOOST_TEST(rbounds->rMax() == rbounds_act->rMax());
    BOOST_TEST(rbounds->averagePhi() == rbounds_act->averagePhi());
    BOOST_TEST(rbounds->halfPhiSector() == rbounds_act->halfPhiSector());

    double rMin = 1, rMax = 5, minHalfX = 2, maxHalfX = 4, stereo = M_PI / 8.;
    auto   dtbounds = std::make_shared<const DiscTrapezoidalBounds>(
        minHalfX, maxHalfX, rMin, rMax, avgPhi, stereo);

    DiscSurface  dtdisc(trf, dtbounds);
    variant_data var_dtdisc = dtdisc.toVariantData();
    std::cout << var_dtdisc;

    DiscSurface dtdisc2(var_dtdisc);
    ;

    BOOST_TEST(dtdisc.transform().isApprox(dtdisc2.transform()));
    BOOST_TEST(dtdisc2.transform().isApprox(*trf));

    const DiscTrapezoidalBounds* dtbounds_act
        = dynamic_cast<const DiscTrapezoidalBounds*>(&dtdisc2.bounds());
    BOOST_TEST(dtbounds->rMin() == dtbounds_act->rMin());
    BOOST_TEST(dtbounds->rMax() == dtbounds_act->rMax());
    BOOST_TEST(dtbounds->minHalflengthX() == dtbounds_act->minHalflengthX());
    BOOST_TEST(dtbounds->maxHalflengthX() == dtbounds_act->maxHalflengthX());
    BOOST_TEST(dtbounds->stereo() == dtbounds_act->stereo());
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
