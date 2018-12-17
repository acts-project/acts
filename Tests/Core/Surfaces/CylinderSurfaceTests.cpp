// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE CylinderSurface Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include <limits>
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
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
    BOOST_TEST(
        Surface::makeShared<CylinderSurface>(pNullTransform, radius, halfZ)
            ->type()
        == Surface::Cylinder);
    BOOST_TEST(
        Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ)->type()
        == Surface::Cylinder);
    //
    /// Constructor with transform pointer, radius, halfZ and halfPhiSector
    BOOST_TEST(Surface::makeShared<CylinderSurface>(
                   pTransform, radius, halfPhiSector, halfZ)
                   ->type()
               == Surface::Cylinder);

    /// Constructor with transform and CylinderBounds pointer
    auto pCylinderBounds
        = std::make_shared<const CylinderBounds>(radius, halfZ);
    BOOST_TEST(Surface::makeShared<CylinderSurface>(pTransform, pCylinderBounds)
                   ->type()
               == Surface::Cylinder);
    //
    //
    /// Copy constructor
    auto cylinderSurfaceObject
        = Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);
    auto copiedCylinderSurface
        = Surface::makeShared<CylinderSurface>(*cylinderSurfaceObject);
    BOOST_TEST(copiedCylinderSurface->type() == Surface::Cylinder);
    BOOST_TEST(*copiedCylinderSurface == *cylinderSurfaceObject);
    //
    /// Copied and transformed
    auto copiedTransformedCylinderSurface
        = Surface::makeShared<CylinderSurface>(*cylinderSurfaceObject,
                                               *pTransform);
    BOOST_TEST(copiedTransformedCylinderSurface->type() == Surface::Cylinder);

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
    BOOST_TEST(pClonedCylinderSurface->type() == Surface::Cylinder);
    //
    /// Test type (redundant)
    BOOST_TEST(cylinderSurfaceObject->type() == Surface::Cylinder);
    //
    /// Test binningPosition
    Vector3D binningPosition{0., 1., 2.};
    BOOST_TEST(cylinderSurfaceObject->binningPosition(BinningValue::binPhi)
               == binningPosition);
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
    BOOST_TEST(cylinderSurfaceObject->referenceFrame(globalPosition, momentum)
                   .isApprox(expectedFrame));
    // check with shift and different momentum
    BOOST_TEST(cylinderSurfaceObject->referenceFrame(globalPositionZ, momentum2)
                   .isApprox(expectedFrame));
    //
    /// Test normal, given 3D position
    Vector3D origin{0., 0., 0.};
    Vector3D normal3D = {0., -1., 0.};
    BOOST_TEST(cylinderSurfaceObject->normal(origin) == normal3D);

    Vector3D pos45deg    = {rootHalf, 1 + rootHalf, 0.};
    Vector3D pos45degZ   = {rootHalf, 1 + rootHalf, 4.};
    Vector3D normal45deg = {rootHalf, rootHalf, 0.};
    // test the normal vector
    BOOST_TEST(cylinderSurfaceObject->normal(pos45deg).isApprox(normal45deg)
               == true);
    // thest that the normal vector is independent of z coordinate
    BOOST_TEST(cylinderSurfaceObject->normal(pos45degZ).isApprox(normal45deg)
               == true);
    //
    /// Test normal given 2D rphi position
    Vector2D positionPiBy2(1.0, 0.);
    Vector3D normalAtPiBy2{std::cos(1.), std::sin(1.), 0.};
    BOOST_TEST(cylinderSurfaceObject->normal(positionPiBy2) == normalAtPiBy2);

    //
    /// Test rotational symmetry axis
    Vector3D symmetryAxis{0., 0., 1.};
    BOOST_TEST(cylinderSurfaceObject->rotSymmetryAxis() == symmetryAxis);
    //
    /// Test bounds
    BOOST_TEST(cylinderSurfaceObject->bounds().type()
               == SurfaceBounds::Cylinder);
    //
    /// Test localToGlobal
    Vector2D localPosition{0., 0.};
    cylinderSurfaceObject->localToGlobal(
        localPosition, momentum, globalPosition);
    Vector3D expectedPosition{1, 1, 2};
    BOOST_TEST(globalPosition == expectedPosition, "Testing localToGlobal");
    //
    /// Testing globalToLocal
    cylinderSurfaceObject->globalToLocal(
        globalPosition, momentum, localPosition);
    Vector2D expectedLocalPosition{0., 0.};
    BOOST_TEST(localPosition == expectedLocalPosition, "Testing globalToLocal");
    //
    /// Test isOnSurface
    Vector3D offSurface{100, 1, 2};
    BOOST_TEST(
        cylinderSurfaceObject->isOnSurface(globalPosition, momentum, true));
    BOOST_TEST(cylinderSurfaceObject->isOnSurface(offSurface, momentum, true)
               == false);
    //
    /// intersectionEstimate
    Vector3D direction{-1., 0, 0};
    auto     intersect
        = cylinderSurfaceObject->intersectionEstimate(offSurface, direction);
    Intersection expectedIntersect{Vector3D{1, 1, 2}, 99., true, 0};
    BOOST_TEST(intersect.valid);
    BOOST_TEST(intersect.position.isApprox(expectedIntersect.position));
    BOOST_TEST(intersect.pathLength == expectedIntersect.pathLength);
    BOOST_TEST(intersect.distance == expectedIntersect.distance);
    //
    /// Test pathCorrection
    BOOST_TEST(cylinderSurfaceObject->pathCorrection(offSurface, momentum)
                   == std::sqrt(3.),
               tt::tolerance(0.01));
    //
    /// Test name
    BOOST_TEST(cylinderSurfaceObject->name()
               == std::string("Acts::CylinderSurface"));
    //
    /// Test dump
    boost::test_tools::output_test_stream dumpOuput;
    cylinderSurfaceObject->dump(dumpOuput);
    BOOST_TEST(dumpOuput.is_equal("Acts::CylinderSurface\n\
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
    BOOST_TEST(*cylinderSurfaceObject == *cylinderSurfaceObject2);
    //
    BOOST_TEST_CHECKPOINT(
        "Create and then assign a CylinderSurface object to the existing one");
    /// Test assignment
    auto assignedCylinderSurface
        = Surface::makeShared<CylinderSurface>(nullptr, 6.6, 5.4);
    *assignedCylinderSurface = *cylinderSurfaceObject;
    /// Test equality of assigned to original
    BOOST_TEST(*assignedCylinderSurface == *cylinderSurfaceObject);
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
    BOOST_TEST(bounds_pl.get<double>("radius") == radius);
    BOOST_TEST(bounds_pl.get<double>("avgPhi") == 0);
    BOOST_TEST(bounds_pl.get<double>("halfPhi") == M_PI);
    BOOST_TEST(bounds_pl.get<double>("halfZ") == halfZ);

    auto cylinder2 = Surface::makeShared<CylinderSurface>(var_cyl);
    auto cylbounds = dynamic_cast<const CylinderBounds*>(&cylinder2->bounds());
    BOOST_TEST(cylbounds->r() == radius);
    BOOST_TEST(cylbounds->halflengthZ() == halfZ);
    BOOST_TEST(cylbounds->halfPhiSector() == M_PI);
    BOOST_TEST(cylbounds->averagePhi() == 0);
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
