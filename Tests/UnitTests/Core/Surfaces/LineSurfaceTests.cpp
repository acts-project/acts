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
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/LineSurfaceStub.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace utf = boost::unit_test;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

// using boost::test_tools::output_test_stream;

BOOST_AUTO_TEST_SUITE(Surfaces)
/// Unit test for creating compliant/non-compliant LineSurface object
BOOST_AUTO_TEST_CASE(LineSurface_Constructors_test) {
  // Default ctor is deleted
  // LineSurfaceStub l;
  // ctor with translation, radius, halfz
  Translation3D translation{0., 1., 2.};
  Transform3D transform(translation);
  auto pTransform = std::make_shared<const Transform3D>(translation);
  const double radius{2.0}, halfz{20.};
  BOOST_CHECK(LineSurfaceStub(pTransform, radius, halfz).constructedOk());
  // ctor with nullptr for LineBounds
  BOOST_CHECK(LineSurfaceStub(pTransform).constructedOk());
  // ctor with LineBounds
  auto pLineBounds = std::make_shared<const LineBounds>(2., 10.0);
  BOOST_CHECK(LineSurfaceStub(pTransform, pLineBounds).constructedOk());
  // ctor with LineBounds, detector element, Identifier
  MaterialProperties properties{1., 1., 1., 20., 10, 5.};
  auto pMaterial =
      std::make_shared<const HomogeneousSurfaceMaterial>(properties);
  DetectorElementStub detElement{pTransform, pLineBounds, 0.2, pMaterial};
  BOOST_CHECK(LineSurfaceStub(pLineBounds, detElement).constructedOk());
  LineSurfaceStub lineToCopy(pTransform, 2.0, 20.);
  // Copy ctor
  BOOST_CHECK(LineSurfaceStub(lineToCopy).constructedOk());
  // Copied and transformed ctor
  BOOST_CHECK(
      LineSurfaceStub(tgContext, lineToCopy, transform).constructedOk());

  /// Construct with nullptr bounds
  DetectorElementStub detElem;
  BOOST_CHECK_THROW(LineSurfaceStub nullBounds(nullptr, detElem),
                    AssertionFailureException);

  BOOST_TEST_MESSAGE(
      "All LineSurface constructors are callable without problem");
}
/// Unit tests of all named methods
BOOST_AUTO_TEST_CASE(LineSurface_allNamedMethods_test) {
  // binningPosition()
  Translation3D translation{0., 1., 2.};
  Transform3D transform(translation);
  auto pTransform = std::make_shared<const Transform3D>(translation);
  LineSurfaceStub line(pTransform, 2.0, 20.);
  Vector3D referencePosition{0., 1., 2.};
  CHECK_CLOSE_ABS(referencePosition, line.binningPosition(tgContext, binX),
                  1e-6);
  //
  // bounds()
  auto pLineBounds = std::make_shared<const LineBounds>(2., 10.0);
  LineSurfaceStub boundedLine(pTransform, pLineBounds);
  const LineBounds& bounds =
      dynamic_cast<const LineBounds&>(boundedLine.bounds());
  BOOST_CHECK_EQUAL(bounds, LineBounds(2., 10.0));
  //
  // globalToLocal()
  Vector3D gpos{0., 1., 0.};
  const Vector3D mom{20., 0., 0.};  // needs more realistic parameters
  Vector2D localPosition;
  BOOST_CHECK(line.globalToLocal(tgContext, gpos, mom, localPosition));
  const Vector2D expectedResult{0, -2};
  CHECK_CLOSE_ABS(expectedResult, localPosition, 1e-6);
  //
  // intersectionEstimate
  const Vector3D direction{0., 1., 2.};
  BoundaryCheck bcheck(false);
  auto intersection = line.intersectionEstimate(tgContext, {0., 0., 0.},
                                                direction.normalized(), bcheck);
  BOOST_CHECK(bool(intersection));
  Vector3D expectedIntersection(0, 1., 2.);
  CHECK_CLOSE_ABS(intersection.position, expectedIntersection,
                  1e-6);  // need more tests..
  //
  // isOnSurface
  const Vector3D insidePosition{0., 2.5, 0.};
  BOOST_CHECK(line.isOnSurface(tgContext, insidePosition, mom,
                               false));  // need better test here
  const Vector3D outsidePosition{100., 100., 200.};
  BOOST_CHECK(!line.isOnSurface(tgContext, outsidePosition, mom, true));
  //
  // localToGlobal
  Vector3D returnedGlobalPosition{0., 0., 0.};
  // Vector2D localPosition{0., 0.};
  const Vector3D momentum{300., 200., 0.};  // find better values!
  line.localToGlobal(tgContext, localPosition, momentum,
                     returnedGlobalPosition);
  const Vector3D expectedGlobalPosition{0, 1, 0};
  CHECK_CLOSE_ABS(returnedGlobalPosition, expectedGlobalPosition, 1e-6);
  //
  // referenceFrame
  Vector3D globalPosition{0., 0., 0.};
  auto returnedRotationMatrix =
      line.referenceFrame(tgContext, globalPosition, momentum);
  double v0 = std::cos(std::atan(2. / 3.));
  double v1 = std::sin(std::atan(2. / 3.));
  RotationMatrix3D expectedRotationMatrix;
  expectedRotationMatrix << -v1, 0., v0, v0, 0., v1, 0., 1., -0.;
  // std::cout<<returnedRotationMatrix<<std::endl;
  // std::cout<<expectedRotationMatrix<<std::endl;
  CHECK_CLOSE_OR_SMALL(returnedRotationMatrix, expectedRotationMatrix, 1e-6,
                       1e-9);
  //
  // name()
  boost::test_tools::output_test_stream output;
  output << line.name();
  BOOST_CHECK(output.is_equal("Acts::LineSurface"));
  //
  // normal
  Vector3D normalVector{0., 0., 1.};  // line direction is same as normal????
  CHECK_CLOSE_ABS(line.normal(tgContext), normalVector, 1e-6);
  //
  // pathCorrection
  auto any3DVector = normalVector;
  CHECK_CLOSE_REL(line.pathCorrection(tgContext, any3DVector, any3DVector), 1.,
                  1e-6);
}
/// Unit test for testing LineSurface assignment
BOOST_AUTO_TEST_CASE(LineSurface_assignment_test) {
  Translation3D translation{0., 1., 2.};
  Transform3D transform(translation);
  auto pTransform = std::make_shared<const Transform3D>(translation);
  LineSurfaceStub originalLine(pTransform, 2.0, 20.);
  LineSurfaceStub assignedLine(pTransform, 1.0, 1.0);
  BOOST_CHECK(assignedLine != originalLine);  // operator != from base
  assignedLine = originalLine;
  BOOST_CHECK(assignedLine == originalLine);  // operator == from base
}

/// Unit test for testing LineSurface alignment derivatives
BOOST_AUTO_TEST_CASE(LineSurfaceAlignment) {
  Translation3D translation{0., 1., 2.};
  Transform3D transform(translation);
  auto pTransform = std::make_shared<const Transform3D>(translation);
  LineSurfaceStub line(pTransform, 2.0, 20.);

  const auto& rotation = pTransform->rotation();
  // The local frame z axis
  const Vector3D localZAxis = rotation.col(2);
  // Check the local z axis is aligned to global z axis
  CHECK_CLOSE_ABS(localZAxis, Vector3D(0., 0., 1.), 1e-15);

  /// Define the track (global) position and direction
  Vector3D globalPosition{1, 2, 4};
  Vector3D momentum{-1, 1, 1};
  Vector3D direction = momentum.normalized();

  // Call the function to calculate the derivative of local frame axes w.r.t its
  // rotation
  const auto& [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);

  // (a) Test the derivative of path length w.r.t. alignment parameters
  const AlignmentRowVector& alignToPath = line.alignmentToPathDerivative(
      tgContext, rotToLocalZAxis, globalPosition, direction);
  // The expected results
  AlignmentRowVector expAlignToPath = AlignmentRowVector::Zero();
  const double value = std::sqrt(3) / 2;
  expAlignToPath << -value, value, 0, -3 * value, -value, 0;
  // Check if the calculated derivative is as expected
  CHECK_CLOSE_ABS(alignToPath, expAlignToPath, 1e-10);

  // (b) Test the derivative of bound track parameters local position w.r.t.
  // position in local 3D Cartesian coordinates
  const auto& loc3DToLocBound =
      line.localCartesianToBoundLocalDerivative(tgContext, globalPosition);
  // Check if the result is as expected
  LocalCartesianToBoundLocalMatrix expLoc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Zero();
  expLoc3DToLocBound << 1 / std::sqrt(2), 1 / std::sqrt(2), 0, 0, 0, 1;
  CHECK_CLOSE_ABS(loc3DToLocBound, expLoc3DToLocBound, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
