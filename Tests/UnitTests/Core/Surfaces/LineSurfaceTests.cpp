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

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/LineSurfaceStub.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

namespace utf = boost::unit_test;

namespace Acts {
class AssertionFailureException;

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
  Translation3 translation{0., 1., 2.};
  Transform3 transform(translation);
  auto pTransform = Transform3(translation);
  const double radius{2.0}, halfz{20.};
  BOOST_CHECK(LineSurfaceStub(pTransform, radius, halfz).constructedOk());
  // ctor with nullptr for LineBounds
  BOOST_CHECK(LineSurfaceStub(pTransform).constructedOk());
  // ctor with LineBounds
  auto pLineBounds = std::make_shared<const LineBounds>(2., 10.0);
  BOOST_CHECK(LineSurfaceStub(pTransform, pLineBounds).constructedOk());
  // ctor with LineBounds, detector element, Identifier
  auto pMaterial =
      std::make_shared<const HomogeneousSurfaceMaterial>(makePercentSlab());
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
  Translation3 translation{0., 1., 2.};
  Transform3 transform(translation);
  LineSurfaceStub line(transform, 2.0, 20.);
  Vector3 referencePosition{0., 1., 2.};
  CHECK_CLOSE_ABS(referencePosition, line.binningPosition(tgContext, binX),
                  1e-6);
  //
  // bounds()
  auto pLineBounds = std::make_shared<const LineBounds>(2., 10.0);
  LineSurfaceStub boundedLine(transform, pLineBounds);
  const LineBounds& bounds =
      dynamic_cast<const LineBounds&>(boundedLine.bounds());
  BOOST_CHECK_EQUAL(bounds, LineBounds(2., 10.0));
  //
  // globalToLocal()
  Vector3 gpos{0., 1., 0.};
  const Vector3 mom{20., 0., 0.};  // needs more realistic parameters
  Vector2 localPosition = line.globalToLocal(tgContext, gpos, mom).value();
  const Vector2 expectedResult{0, -2};
  CHECK_CLOSE_ABS(expectedResult, localPosition, 1e-6);
  //
  // intersection
  const Vector3 direction{0., 1., 2.};
  BoundaryCheck bcheck(false);
  auto sfIntersection =
      line.intersect(tgContext, {0., 0., 0.}, direction.normalized(), bcheck);
  BOOST_CHECK(bool(sfIntersection));
  Vector3 expectedIntersection(0, 1., 2.);
  CHECK_CLOSE_ABS(sfIntersection.intersection.position, expectedIntersection,
                  1e-6);  // need more tests..
  BOOST_CHECK_EQUAL(sfIntersection.object, &line);
  //
  // isOnSurface
  const Vector3 insidePosition{0., 2.5, 0.};
  BOOST_CHECK(line.isOnSurface(tgContext, insidePosition, mom,
                               false));  // need better test here
  const Vector3 outsidePosition{100., 100., 200.};
  BOOST_CHECK(!line.isOnSurface(tgContext, outsidePosition, mom, true));
  //
  // localToGlobal
  Vector3 returnedGlobalPosition{0., 0., 0.};
  // Vector2 localPosition{0., 0.};
  const Vector3 momentum{300., 200., 0.};  // find better values!
  returnedGlobalPosition =
      line.localToGlobal(tgContext, localPosition, momentum);
  const Vector3 expectedGlobalPosition{0, 1, 0};
  CHECK_CLOSE_ABS(returnedGlobalPosition, expectedGlobalPosition, 1e-6);
  //
  // referenceFrame
  Vector3 globalPosition{0., 0., 0.};
  auto returnedRotationMatrix =
      line.referenceFrame(tgContext, globalPosition, momentum);
  double v0 = std::cos(std::atan(2. / 3.));
  double v1 = std::sin(std::atan(2. / 3.));
  RotationMatrix3 expectedRotationMatrix;
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
  Vector3 normalVector{0., 0., 1.};  // line direction is same as normal????
  CHECK_CLOSE_ABS(line.normal(tgContext), normalVector, 1e-6);
  //
  // pathCorrection
  auto any3DVector = normalVector;
  CHECK_CLOSE_REL(line.pathCorrection(tgContext, any3DVector, any3DVector), 1.,
                  1e-6);
}
/// Unit test for testing LineSurface assignment
BOOST_AUTO_TEST_CASE(LineSurface_assignment_test) {
  Translation3 translation{0., 1., 2.};
  Transform3 transform(translation);
  LineSurfaceStub originalLine(transform, 2.0, 20.);
  LineSurfaceStub assignedLine(transform, 1.0, 1.0);
  BOOST_CHECK(assignedLine != originalLine);  // operator != from base
  assignedLine = originalLine;
  BOOST_CHECK(assignedLine == originalLine);  // operator == from base
}

/// Unit test for testing LineSurface alignment derivatives
BOOST_AUTO_TEST_CASE(LineSurfaceAlignment) {
  Translation3 translation{0., 1., 2.};
  Transform3 transform(translation);
  LineSurfaceStub line(transform, 2.0, 20.);

  const auto& rotation = transform.rotation();
  // The local frame z axis
  const Vector3 localZAxis = rotation.col(2);
  // Check the local z axis is aligned to global z axis
  CHECK_CLOSE_ABS(localZAxis, Vector3(0., 0., 1.), 1e-15);

  // Define the track (global) position and direction
  Vector3 globalPosition{1, 2, 4};
  Vector3 momentum{-1, 1, 1};
  Vector3 direction = momentum.normalized();
  // Construct a free parameters
  FreeVector parameters = FreeVector::Zero();
  parameters.head<3>() = globalPosition;
  parameters.segment<3>(eFreeDir0) = direction;

  // (a) Test the derivative of path length w.r.t. alignment parameters
  const AlignmentToPathMatrix& alignToPath =
      line.alignmentToPathDerivative(tgContext, parameters);
  // The expected results
  AlignmentToPathMatrix expAlignToPath = AlignmentToPathMatrix::Zero();
  const double value = std::sqrt(3) / 2;
  expAlignToPath << -value, value, 0, -3 * value, -value, 0;
  // Check if the calculated derivative is as expected
  CHECK_CLOSE_ABS(alignToPath, expAlignToPath, 1e-10);

  // (b) Test the derivative of bound track parameters local position w.r.t.
  // position in local 3D Cartesian coordinates
  const auto& loc3DToLocBound =
      line.localCartesianToBoundLocalDerivative(tgContext, globalPosition);
  // Check if the result is as expected
  ActsMatrix<2, 3> expLoc3DToLocBound = ActsMatrix<2, 3>::Zero();
  expLoc3DToLocBound << 1 / std::sqrt(2), 1 / std::sqrt(2), 0, 0, 0, 1;
  CHECK_CLOSE_ABS(loc3DToLocBound, expLoc3DToLocBound, 1e-10);
}

BOOST_AUTO_TEST_CASE(LineSurfaceTransformRoundTrip) {
  LineSurfaceStub surface(Transform3::Identity());

  auto roundTrip = [&surface](const Vector3& pos, const Vector3& dir) {
    auto intersection = surface.intersect(tgContext, pos, dir);
    Vector3 global = intersection.intersection.position;
    Vector2 local = *surface.globalToLocal(tgContext, global, dir);
    Vector3 global2 = surface.localToGlobal(tgContext, local, dir);
    return std::make_tuple(global, local, global2);
  };

  {
    Vector3 pos = {-0.02801, 0.00475611, 0.285106};
    Vector3 dir = Vector3(-0.03951, -0.221457, -0.564298).normalized();

    auto [global, local, global2] = roundTrip(pos, dir);

    CHECK_CLOSE_ABS(global, global2, 1e-10);
  }

  {
    Vector3 pos = {-64.2892, 65.2697, -0.839014};
    Vector3 dir = Vector3(-0.236602, -0.157616, 0.956786).normalized();

    auto [global, local, global2] = roundTrip(pos, dir);

    CHECK_CLOSE_ABS(global, global2, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(LineSurfaceTransformRoundTripEtaStability) {
  LineSurfaceStub surface(Transform3::Identity());

  // eta=6 is already crashing
  const std::vector<double> etas = {0, 1, 2, 3, 4, 5};

  for (double eta : etas) {
    Vector3 pca = {5, 0, 0};
    Vector3 dir = makeDirectionUnitFromPhiEta(M_PI_2, eta);
    Vector3 pos = pca + dir;

    auto intersection = surface.intersect(tgContext, pos, dir);

    Vector3 global = intersection.intersection.position;
    Vector2 local = *surface.globalToLocal(tgContext, global, dir);
    Vector3 global2 = surface.localToGlobal(tgContext, local, dir);

    CHECK_CLOSE_ABS(global, global2, 1e-10);
    CHECK_CLOSE_ABS(pca, global2, 1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
