// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
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
#include "ACTS/Material/HomogeneousSurfaceMaterial.hpp"
#include "ACTS/Surfaces/LineSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
//
#include "DetectorElementStub.hpp"
#include "LineSurfaceStub.hpp"
//
#include <limits>

namespace utf    = boost::unit_test;
const double inf = std::numeric_limits<double>::infinity();
const double NaN = std::numeric_limits<double>::quiet_NaN();

namespace Acts {

namespace Test {
  // using boost::test_tools::output_test_stream;

  BOOST_AUTO_TEST_SUITE(Surfaces);
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
    auto pLineBounds = std::make_shared<LineBounds>(10.0);
    BOOST_CHECK(LineSurfaceStub(pTransform, pLineBounds).constructedOk());
    // ctor with LineBounds, detector element, Identifier
    Identifier         identifier{2};
    MaterialProperties properties{1., 1., 1., 20., 10, 5.};
    auto               pMaterial
        = std::make_shared<const HomogeneousSurfaceMaterial>(properties);
    DetectorElementStub detElement{
        identifier, pTransform, pLineBounds, 0.2, pMaterial};
    BOOST_CHECK(
        LineSurfaceStub(pLineBounds, detElement, identifier).constructedOk());
    LineSurfaceStub lineToCopy(pTransform, 2.0, 20.);
    // Copy ctor
    BOOST_CHECK(LineSurfaceStub(lineToCopy).constructedOk());
    // Copied and transformed ctor
    BOOST_CHECK(LineSurfaceStub(lineToCopy, transform).constructedOk());
    BOOST_TEST_MESSAGE(
        "All LineSurface constructors are callable without problem");
  }
  /// Unit tests of all named methods
  BOOST_AUTO_TEST_CASE(LineSurface_allNamedMethods_test)
  {
    // binningPosition()
    Translation3D   translation{0., 1., 2.};
    Transform3D     transform(translation);
    auto pTransform = std::make_shared<const Transform3D>(translation);
    LineSurfaceStub line(pTransform, 2.0, 20.);
    Vector3D        referencePosition{0., 1., 2.};
    BOOST_TEST(referencePosition == line.binningPosition(binX));
    //
    // bounds()
    auto              pLineBounds = std::make_shared<LineBounds>(10.0);
    LineSurfaceStub   boundedLine(pTransform, pLineBounds);
    const LineBounds& bounds
        = dynamic_cast<const LineBounds&>(boundedLine.bounds());
    BOOST_TEST(bounds == LineBounds(10.0),
               "bounds() equals construction bounds");
    //
    // globalToLocal()
    Vector3D gpos{
        0., 1., 0,
    };
    const Vector3D mom{20., 0., 0.};  // needs more realistic parameters
    Vector2D       localPosition;
    BOOST_CHECK(line.globalToLocal(gpos, mom, localPosition));
    const Vector2D expectedResult{0, -2};
    BOOST_CHECK(expectedResult == localPosition);
    //
    // intersectionEstimate
    const Vector3D direction{0., 1., 2.};
    bool           forceDir = false;
    BoundaryCheck  bcheck(false);
    auto           intersection
        = line.intersectionEstimate({0., 0., 0.}, direction, forceDir, bcheck);
    BOOST_CHECK(intersection.valid);
    Vector3D expectedIntersection(0, -1. / 3., -2. / 3.);
    BOOST_TEST(intersection.position
               == expectedIntersection);  // need more tests..
    //
    // isOnSurface
    const Vector3D insidePosition{0., 2.5, 0.};
    BOOST_TEST(
        line.isOnSurface(insidePosition, false));  // need better test here
    const Vector3D outsidePosition{100., 100., 200.};
    BOOST_TEST(!line.isOnSurface(outsidePosition, true));
    //
    // lineDirection
    const Vector3D zDirection{0., 0., 1.};
    BOOST_TEST(line.lineDirection() == zDirection);
    //
    // localToGlobal
    Vector3D returnedGlobalPosition{0., 0., 0.};
    // Vector2D localPosition{0., 0.};
    const Vector3D momentum{300., 200., 0.};  // find better values!
    line.localToGlobal(localPosition, momentum, returnedGlobalPosition);
    const Vector3D expectedGlobalPosition{0, 1, 0};
    BOOST_TEST(returnedGlobalPosition == expectedGlobalPosition);
    //
    // measurementFrame
    Vector3D globalPosition{0., 0., 0.};
    auto     returnedRotationMatrix
        = line.measurementFrame(globalPosition, momentum);
    double           v0 = std::cos(std::atan(2. / 3.));
    double           v1 = std::sin(std::atan(2. / 3.));
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
    BOOST_TEST(line.normal() == normalVector);
    //
    // pathCorrection
    auto any3DVector = normalVector;
    BOOST_TEST(line.pathCorrection(any3DVector, any3DVector) == 1.);
  }
  /// Unit test for testing LineSurface assignment
  BOOST_AUTO_TEST_CASE(LineSurface_assignment_test)
  {
    Translation3D   translation{0., 1., 2.};
    Transform3D     transform(translation);
    auto pTransform = std::make_shared<const Transform3D>(translation);
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
  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts
