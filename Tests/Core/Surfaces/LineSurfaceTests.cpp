// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Line Surface Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantData.hpp"
//
#include "../Utilities/TestHelper.hpp"
#include "DetectorElementStub.hpp"
#include "LineSurfaceStub.hpp"
//
#include <limits>

namespace utf = boost::unit_test;

namespace Acts {

namespace Test {
  // using boost::test_tools::output_test_stream;

  BOOST_AUTO_TEST_SUITE(Surfaces)
  /// Unit test for creating compliant/non-compliant LineSurface object
  BOOST_AUTO_TEST_CASE(LineSurface_Constructors_test)
  {
    // Default ctor is deleted
    // LineSurfaceStub l;
    // ctor with translation, radius, halfz
    Translation3D translation{0., 1., 2.};
    Transform3D   transform(translation);
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    const double  radius{2.0}, halfz{20.};
    BOOST_CHECK(LineSurfaceStub(pTransform, radius, halfz).constructedOk());
    // ctor with nullptr for LineBounds
    BOOST_CHECK(LineSurfaceStub(pTransform).constructedOk());
    // ctor with LineBounds
    auto pLineBounds = std::make_shared<const LineBounds>(10.0);
    BOOST_CHECK(LineSurfaceStub(pTransform, pLineBounds).constructedOk());
    // ctor with LineBounds, detector element, Identifier
    Identifier         identifier{2};
    MaterialProperties properties{1., 1., 1., 20., 10, 5.};
    auto               pMaterial
        = std::make_shared<const HomogeneousSurfaceMaterial>(properties);
    DetectorElementStub detElement{
        identifier, pTransform, pLineBounds, 0.2, pMaterial};
    BOOST_CHECK(LineSurfaceStub(pLineBounds, detElement).constructedOk());
    LineSurfaceStub lineToCopy(pTransform, 2.0, 20.);
    // Copy ctor
    BOOST_CHECK(LineSurfaceStub(lineToCopy).constructedOk());
    // Copied and transformed ctor
    BOOST_CHECK(LineSurfaceStub(lineToCopy, transform).constructedOk());

    /// Construct with nullptr bounds
    DetectorElementStub detElem;
    BOOST_CHECK_THROW(LineSurfaceStub nullBounds(nullptr, detElem),
                      AssertionFailureException);

    BOOST_TEST_MESSAGE(
        "All LineSurface constructors are callable without problem");
  }
  /// Unit tests of all named methods
  BOOST_AUTO_TEST_CASE(LineSurface_allNamedMethods_test)
  {
    // binningPosition()
    Translation3D translation{0., 1., 2.};
    Transform3D   transform(translation);
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    LineSurfaceStub line(pTransform, 2.0, 20.);
    Vector3D        referencePosition{0., 1., 2.};
    checkCloseVec3D(referencePosition, line.binningPosition(binX));
    //
    // bounds()
    auto              pLineBounds = std::make_shared<const LineBounds>(10.0);
    LineSurfaceStub   boundedLine(pTransform, pLineBounds);
    const LineBounds& bounds
        = dynamic_cast<const LineBounds&>(boundedLine.bounds());
    BOOST_TEST(bounds == LineBounds(10.0),
               "bounds() equals construction bounds");
    //
    // globalToLocal()
    Vector3D       gpos{0., 1., 0.};
    const Vector3D mom{20., 0., 0.};  // needs more realistic parameters
    Vector2D       localPosition;
    BOOST_CHECK(line.globalToLocal(gpos, mom, localPosition));
    const Vector2D expectedResult{0, -2};
    checkCloseVec2D(expectedResult, localPosition);
    //
    // intersectionEstimate
    const Vector3D      direction{0., 1., 2.};
    NavigationDirection navDir = anyDirection;
    BoundaryCheck       bcheck(false);
    auto                intersection = line.intersectionEstimate(
        {0., 0., 0.}, direction.normalized(), navDir, bcheck);
    BOOST_CHECK(intersection.valid);
    Vector3D expectedIntersection(0, 1., 2.);
    checkCloseVec3D(intersection.position,
                    expectedIntersection);  // need more tests..
    //
    // isOnSurface
    const Vector3D insidePosition{0., 2.5, 0.};
    BOOST_TEST(
        line.isOnSurface(insidePosition, mom, false));  // need better test here
    const Vector3D outsidePosition{100., 100., 200.};
    BOOST_TEST(!line.isOnSurface(outsidePosition, mom, true));
    //
    // lineDirection
    const Vector3D zDirection{0., 0., 1.};
    checkCloseVec3D(line.lineDirection(), zDirection);
    //
    // localToGlobal
    Vector3D returnedGlobalPosition{0., 0., 0.};
    // Vector2D localPosition{0., 0.};
    const Vector3D momentum{300., 200., 0.};  // find better values!
    line.localToGlobal(localPosition, momentum, returnedGlobalPosition);
    const Vector3D expectedGlobalPosition{0, 1, 0};
    checkCloseVec3D(returnedGlobalPosition, expectedGlobalPosition);
    //
    // referenceFrame
    Vector3D globalPosition{0., 0., 0.};
    auto returnedRotationMatrix = line.referenceFrame(globalPosition, momentum);
    double           v0         = std::cos(std::atan(2. / 3.));
    double           v1         = std::sin(std::atan(2. / 3.));
    RotationMatrix3D expectedRotationMatrix;
    expectedRotationMatrix << -v1, 0., v0, v0, 0., v1, 0., 1., -0.;
    // std::cout<<returnedRotationMatrix<<std::endl;
    // std::cout<<expectedRotationMatrix<<std::endl;
    BOOST_TEST(returnedRotationMatrix.isApprox(expectedRotationMatrix));
    //
    // name()
    boost::test_tools::output_test_stream output;
    output << line.name();
    BOOST_TEST(output.is_equal("Acts::LineSurface"));
    //
    // normal
    Vector3D normalVector{0., 0., 1.};  // line direction is same as normal????
    checkCloseVec3D(line.normal(), normalVector);
    //
    // pathCorrection
    auto any3DVector = normalVector;
    BOOST_CHECK_CLOSE_FRACTION(
        line.pathCorrection(any3DVector, any3DVector), 1., 1e-6);
  }
  /// Unit test for testing LineSurface assignment
  BOOST_AUTO_TEST_CASE(LineSurface_assignment_test)
  {
    Translation3D translation{0., 1., 2.};
    Transform3D   transform(translation);
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    LineSurfaceStub originalLine(pTransform, 2.0, 20.);
    LineSurfaceStub assignedLine(pTransform, 1.0, 1.0);
    BOOST_TEST(assignedLine != originalLine,
               "LineSurfaces are different before assignment");  // operator !=
                                                                 // from base
    assignedLine = originalLine;
    BOOST_TEST(assignedLine == originalLine,
               "LineSurfaces are equal value after assignment");  // operator ==
                                                                  // from base
  }

  BOOST_AUTO_TEST_CASE(LineSurface_toVariantData)
  {
    double        radius = 2.0, hlZ = 20;
    Translation3D translation{0., 1., 2.};
    Transform3D   transform(translation);
    auto          pTransform = std::make_shared<const Transform3D>(translation);
    LineSurfaceStub line(pTransform, radius, hlZ);
    variant_data    var_line = line.toVariantData();
    std::cout << var_line << std::endl;

    const variant_map& pl
        = boost::get<variant_map>(var_line).get<variant_map>("payload");
    const variant_map& bounds_pl
        = pl.get<variant_map>("bounds").get<variant_map>("payload");
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("radius"), radius);
    BOOST_CHECK_EQUAL(bounds_pl.get<double>("halfZ"), hlZ);

    LineSurfaceStub line2(var_line);
    auto            lbounds = dynamic_cast<const LineBounds*>(&line2.bounds());
    BOOST_CHECK_EQUAL(lbounds->r(), radius);
    BOOST_CHECK_EQUAL(lbounds->halflengthZ(), hlZ);
  }

  BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
