// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE CylinderSurface Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"

namespace tt = boost::test_tools;
using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(CylinderSurfaces)
  /// Unit test for creating compliant/non-compliant CylinderSurface object
  BOOST_AUTO_TEST_CASE(CylinderSurfaceConstruction)
  {
    // CylinderSurface default constructor is deleted
    //
    /// Constructor with transform pointer, null or valid, radius and halfZ
    double        radius(1.0), halfZ(10.), halfPhiSector(M_PI / 8.);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    auto          pNullTransform = std::make_shared<const Transform3D>();
    BOOST_CHECK_EQUAL(
        Surface::makeShared<CylinderSurface>(pNullTransform, radius, halfZ)
            ->type(),
        Surface::Cylinder);
    BOOST_CHECK_EQUAL(
        Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ)->type(),
        Surface::Cylinder);
    //
    /// Constructor with transform pointer, radius, halfZ and halfPhiSector
    BOOST_CHECK_EQUAL(Surface::makeShared<CylinderSurface>(
                          pTransform, radius, halfPhiSector, halfZ)
                          ->type(),
                      Surface::Cylinder);

    /// Constructor with transform and CylinderBounds pointer
    auto pCylinderBounds
        = std::make_shared<const CylinderBounds>(radius, halfZ);
    BOOST_CHECK_EQUAL(
        Surface::makeShared<CylinderSurface>(pTransform, pCylinderBounds)
            ->type(),
        Surface::Cylinder);
    //
    //
    /// Copy constructor
    auto cylinderSurfaceObject
        = Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
    auto copiedCylinderSurface
        = Surface::makeShared<CylinderSurface>(*cylinderSurfaceObject);
    BOOST_CHECK_EQUAL(copiedCylinderSurface->type(), Surface::Cylinder);
    BOOST_CHECK_EQUAL(*copiedCylinderSurface, *cylinderSurfaceObject);
    //
    /// Copied and transformed
    auto copiedTransformedCylinderSurface
        = Surface::makeShared<CylinderSurface>(*cylinderSurfaceObject,
                                               *pTransform);
    BOOST_CHECK_EQUAL(copiedTransformedCylinderSurface->type(),
                      Surface::Cylinder);

    /// Construct with nullptr bounds
    BOOST_CHECK_THROW(auto nullBounds
                      = Surface::makeShared<CylinderSurface>(nullptr, nullptr),
                      AssertionFailureException);
  }
  //
  /// Unit test for testing CylinderSurface properties
  BOOST_AUTO_TEST_CASE(CylinderSurfaceProperties)
  {
    /// Test clone method
    double        radius(1.0), halfZ(10.);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    // auto pNullTransform = std::make_shared<const Transform3D>();
    auto cylinderSurfaceObject
        = Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
    //
    auto pClonedCylinderSurface = cylinderSurfaceObject->clone();
    BOOST_CHECK_EQUAL(pClonedCylinderSurface->type(), Surface::Cylinder);
    //
    /// Test type (redundant)
    BOOST_CHECK_EQUAL(cylinderSurfaceObject->type(), Surface::Cylinder);
    //
    /// Test binningPosition
    Vector3D binningPosition{0., 1., 2.};
    CHECK_CLOSE_ABS(
        cylinderSurfaceObject->binningPosition(BinningValue::binPhi),
        binningPosition,
        1e-9);
    //
    /// Test referenceFrame
    double           rootHalf = std::sqrt(0.5);
    Vector3D         globalPosition{rootHalf, 1. - rootHalf, 0.};
    Vector3D         globalPositionZ{rootHalf, 1. - rootHalf, 2.0};
    Vector3D         momentum{15., 15., 15.};
    Vector3D         momentum2{6.6, -3., 2.};
    RotationMatrix3D expectedFrame;
    expectedFrame << rootHalf, 0., rootHalf, rootHalf, 0., -rootHalf, 0., 1.,
        0.;
    // check without shift
    CHECK_CLOSE_OR_SMALL(
        cylinderSurfaceObject->referenceFrame(globalPosition, momentum),
        expectedFrame,
        1e-6,
        1e-9);
    // check with shift and different momentum
    CHECK_CLOSE_OR_SMALL(
        cylinderSurfaceObject->referenceFrame(globalPositionZ, momentum2),
        expectedFrame,
        1e-6,
        1e-9);
    //
    /// Test normal, given 3D position
    Vector3D origin{0., 0., 0.};
    Vector3D normal3D = {0., -1., 0.};
    CHECK_CLOSE_ABS(cylinderSurfaceObject->normal(origin), normal3D, 1e-9);

    Vector3D pos45deg    = {rootHalf, 1 + rootHalf, 0.};
    Vector3D pos45degZ   = {rootHalf, 1 + rootHalf, 4.};
    Vector3D normal45deg = {rootHalf, rootHalf, 0.};
    // test the normal vector
    CHECK_CLOSE_ABS(
        cylinderSurfaceObject->normal(pos45deg), normal45deg, 1e-6 * rootHalf);
    // thest that the normal vector is independent of z coordinate
    CHECK_CLOSE_ABS(
        cylinderSurfaceObject->normal(pos45degZ), normal45deg, 1e-6 * rootHalf);
    //
    /// Test normal given 2D rphi position
    Vector2D positionPiBy2(1.0, 0.);
    Vector3D normalAtPiBy2{std::cos(1.), std::sin(1.), 0.};
    CHECK_CLOSE_ABS(
        cylinderSurfaceObject->normal(positionPiBy2), normalAtPiBy2, 1e-9);

    //
    /// Test rotational symmetry axis
    Vector3D symmetryAxis{0., 0., 1.};
    CHECK_CLOSE_ABS(
        cylinderSurfaceObject->rotSymmetryAxis(), symmetryAxis, 1e-9);
    //
    /// Test bounds
    BOOST_CHECK_EQUAL(cylinderSurfaceObject->bounds().type(),
                      SurfaceBounds::Cylinder);
    //
    /// Test localToGlobal
    Vector2D localPosition{0., 0.};
    cylinderSurfaceObject->localToGlobal(
        localPosition, momentum, globalPosition);
    Vector3D expectedPosition{1, 1, 2};
    BOOST_CHECK_EQUAL(globalPosition, expectedPosition);
    //
    /// Testing globalToLocal
    cylinderSurfaceObject->globalToLocal(
        globalPosition, momentum, localPosition);
    Vector2D expectedLocalPosition{0., 0.};
    BOOST_CHECK_EQUAL(localPosition, expectedLocalPosition);
    //
    /// Test isOnSurface
    Vector3D offSurface{100, 1, 2};
    BOOST_CHECK(
        cylinderSurfaceObject->isOnSurface(globalPosition, momentum, true));
    BOOST_CHECK(
        !cylinderSurfaceObject->isOnSurface(offSurface, momentum, true));
    //
    /// intersectionEstimate
    Vector3D direction{-1., 0, 0};
    auto     intersect
        = cylinderSurfaceObject->intersectionEstimate(offSurface, direction);
    Intersection expectedIntersect{Vector3D{1, 1, 2}, 99., true, 0};
    BOOST_CHECK(intersect.valid);
    CHECK_CLOSE_ABS(intersect.position, expectedIntersect.position, 1e-9);
    CHECK_CLOSE_ABS(intersect.pathLength, expectedIntersect.pathLength, 1e-9);
    CHECK_CLOSE_ABS(intersect.distance, expectedIntersect.distance, 1e-9);
    //
    /// Test pathCorrection
    CHECK_CLOSE_REL(cylinderSurfaceObject->pathCorrection(offSurface, momentum),
                    std::sqrt(3.),
                    0.01);
    //
    /// Test name
    BOOST_CHECK_EQUAL(cylinderSurfaceObject->name(),
                      std::string("Acts::CylinderSurface"));
    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    cylinderSurfaceObject->dump(dumpOuput);
    BOOST_CHECK(dumpOuput.is_equal("Acts::CylinderSurface\n\
     Center position  (x, y, z) = (0.0000, 1.0000, 2.0000)\n\
     Rotation:             colX = (1.000000, 0.000000, 0.000000)\n\
                           colY = (0.000000, 1.000000, 0.000000)\n\
                           colZ = (0.000000, 0.000000, 1.000000)\n\
     Bounds  : Acts::CylinderBounds: (radius, averagePhi, halfPhiSector, halflengthInZ) = (1.0000000, 0.0000000, 3.1415927, 10.0000000)"));
  }

  BOOST_AUTO_TEST_CASE(EqualityOperators)
  {
    double        radius(1.0), halfZ(10.);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    auto          cylinderSurfaceObject
        = Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
    //
    auto cylinderSurfaceObject2
        = Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
    //
    /// Test equality operator
    BOOST_CHECK_EQUAL(*cylinderSurfaceObject, *cylinderSurfaceObject2);
    //
    BOOST_TEST_CHECKPOINT(
        "Create and then assign a CylinderSurface object to the existing one");
    /// Test assignment
    auto assignedCylinderSurface
        = Surface::makeShared<CylinderSurface>(nullptr, 6.6, 5.4);
    *assignedCylinderSurface = *cylinderSurfaceObject;
    /// Test equality of assigned to original
    BOOST_CHECK_EQUAL(*assignedCylinderSurface, *cylinderSurfaceObject);
  }

  BOOST_AUTO_TEST_CASE(CylinderSurface_toVariantData)
  {
    double        radius(1.0), halfZ(10.);
    Translation3D translation{0., 1., 2.};
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    auto          cylinder
        = Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);

    variant_data var_cyl = cylinder->toVariantData();
    std::cout << var_cyl << std::endl;

    const variant_map& pl
        = boost::get<variant_map>(var_cyl).get<variant_map>("payload");
    const variant_map& bounds_pl
        = pl.get<variant_map>("bounds").get<variant_map>("payload");
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("radius"), radius);
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("avgPhi"), 0);
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("halfPhi"), M_PI);
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("halfZ"), halfZ);

    auto cylinder2 = Surface::makeShared<CylinderSurface>(var_cyl);
    auto cylbounds = dynamic_cast<const CylinderBounds*>(&cylinder2->bounds());
    BOOST_CHECK_EQUAL(cylbounds->r(), radius);
    BOOST_CHECK_EQUAL(cylbounds->halflengthZ(), halfZ);
    BOOST_CHECK_EQUAL(cylbounds->halfPhiSector(), M_PI);
    BOOST_CHECK_EQUAL(cylbounds->averagePhi(), 0);
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
