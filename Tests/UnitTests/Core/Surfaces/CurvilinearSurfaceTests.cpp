// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

BOOST_AUTO_TEST_CASE(jacobian_test) {
  // (1a) Standard test with curvilinear not glazingly close to z axis
  Vector3 direction = Vector3(7., 8., 9.).normalized();
  CurvilinearSurface surface = CurvilinearSurface(direction);
  FreeToBoundMatrix f2cJacobian = surface.freeToBoundJacobian().value();

  double phi = VectorHelpers::phi(direction);
  double theta = VectorHelpers::theta(direction);
  double sinPhi = std::sin(phi);
  double cosPhi = std::cos(phi);
  double sinTheta = std::sin(theta);
  double cosTheta = std::cos(theta);

  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc0, eFreePos0), -sinPhi, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc0, eFreePos1), cosPhi, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc1, eFreePos0), -cosPhi * cosTheta, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc1, eFreePos1), -sinPhi * cosTheta, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc1, eFreePos2), sinTheta, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundTime, eFreeTime), 1., 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundPhi, eFreeDir0), -sinPhi / sinTheta, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundPhi, eFreeDir1), cosPhi / sinTheta, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundTheta, eFreeDir0), cosPhi * cosTheta, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundTheta, eFreeDir1), sinPhi * cosTheta, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundTheta, eFreeDir2), -sinTheta, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundQOverP, eFreeQOverP), 1., 1e-5);

  // (2a) Standard test with curvilinear not glazingly close to z axis
  direction = Vector3(7., 8., 9.).normalized();
  surface = CurvilinearSurface(direction);
  BoundToFreeMatrix c2fJacobian = surface.boundToFreeJacobian();

  phi = VectorHelpers::phi(direction);
  theta = VectorHelpers::theta(direction);
  sinPhi = std::sin(phi);
  cosPhi = std::cos(phi);
  sinTheta = std::sin(theta);
  cosTheta = std::cos(theta);

  CHECK_CLOSE_REL(c2fJacobian(eFreePos0, eBoundLoc0), -sinPhi, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreePos0, eBoundLoc1), -cosPhi * cosTheta, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreePos1, eBoundLoc0), cosPhi, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreePos1, eBoundLoc1), -sinPhi * cosTheta, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreePos2, eBoundLoc1), sinTheta, 1e-5);
  // Time parameter: stays as is
  CHECK_CLOSE_REL(c2fJacobian(eFreeTime, eBoundTime), 1, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreeDir0, eBoundPhi), -sinTheta * sinPhi, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreeDir0, eBoundTheta), cosTheta * cosPhi, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreeDir1, eBoundPhi), sinTheta * cosPhi, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreeDir1, eBoundTheta), cosTheta * sinPhi, 1e-5);
  CHECK_CLOSE_REL(c2fJacobian(eFreeDir2, eBoundTheta), -sinTheta, 1e-5);
  // Q/P parameter: stays as is
  CHECK_CLOSE_REL(c2fJacobian(eFreeQOverP, eBoundQOverP), 1, 1e-5);
}

BOOST_AUTO_TEST_CASE(SynthesizedPlaneUsesSameJacobians) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  for (const Vector3& direction :
       {Vector3(0.3, 0.4, 0.5).normalized().eval(),
        Vector3(1e-10, 0., 1.).normalized().eval(),
        Vector3(1e-10, 0., -1.).normalized().eval()}) {
    const Vector3 position(1., 2., 3.);
    const CurvilinearSurface curvilinear(position, direction);
    const auto plane = curvilinear.planeSurface();
    BOOST_CHECK(curvilinear.boundParametrization().isGlobalSpherical());
    BOOST_CHECK(plane->boundParametrization(gctx).isGlobalSpherical());
    const auto forward = curvilinear.boundToFreeJacobian();
    const auto inverse = curvilinear.freeToBoundJacobian().value();
    CHECK_CLOSE_ABS(
        forward, plane->boundToFreeJacobian(gctx, position, direction), 1e-14);
    CHECK_CLOSE_ABS(
        inverse, plane->freeToBoundJacobian(gctx, position, direction).value(),
        1e-14);
    CHECK_CLOSE_ABS((inverse * forward).eval(), BoundMatrix::Identity(), 1e-14);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
