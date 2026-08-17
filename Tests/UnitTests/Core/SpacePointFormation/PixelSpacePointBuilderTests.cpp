// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/SpacePointFormation/PixelSpacePointBuilder.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <cmath>
#include <memory>
#include <numbers>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

constexpr double sigmaLoc0 = 15_um;
constexpr double sigmaLoc1 = 40_um;

/// A plane whose local axes are given explicitly, centred at @p center
std::shared_ptr<PlaneSurface> makePlane(const Vector3& center,
                                        const Vector3& loc0,
                                        const Vector3& loc1) {
  Transform3 transform = Transform3::Identity();
  transform.matrix().block<3, 1>(0, 0) = loc0.normalized();
  transform.matrix().block<3, 1>(0, 1) = loc1.normalized();
  transform.matrix().block<3, 1>(0, 2) = loc0.cross(loc1).normalized();
  transform.matrix().block<3, 1>(0, 3) = center;
  return Surface::makeShared<PlaneSurface>(
      transform, std::make_shared<RectangleBounds>(50_mm, 50_mm));
}

SquareMatrix2 localCov() {
  SquareMatrix2 cov = SquareMatrix2::Zero();
  cov(0, 0) = sigmaLoc0 * sigmaLoc0;
  cov(1, 1) = sigmaLoc1 * sigmaLoc1;
  return cov;
}

/// The same variance from a numerically differentiated Jacobian
Vector2 numericVarianceZR(const PlaneSurface& surface, const Vector3& position,
                          const SquareMatrix2& cov) {
  const SquareMatrix3 rot =
      surface.referenceFrame(gctx, position, Vector3::UnitZ());
  constexpr double h = 1e-6;

  SquareMatrix2 jac = SquareMatrix2::Zero();
  for (int i = 0; i < 2; ++i) {
    const Vector3 step = h * rot.col(i);
    const Vector3 plus = position + step;
    const Vector3 minus = position - step;
    jac(0, i) = (plus.z() - minus.z()) / (2 * h);
    jac(1, i) =
        (fastHypot(plus.x(), plus.y()) - fastHypot(minus.x(), minus.y())) /
        (2 * h);
  }
  return (jac * cov * jac.transpose()).diagonal();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PixelSpacePointBuilderSuite)

/// A barrel module measures r-phi and z, so all of the variance lands in z
BOOST_AUTO_TEST_CASE(BarrelModule) {
  const double phi = 0.3;
  const double r = 120_mm;
  const Vector3 position(r * std::cos(phi), r * std::sin(phi), 250_mm);
  const Vector3 rPhi(-std::sin(phi), std::cos(phi), 0);
  const Vector3 z = Vector3::UnitZ();

  auto surface = makePlane(position, rPhi, z);
  const Vector2 varZR = PixelSpacePointBuilder::computeVarianceZR(
      gctx, *surface, position, localCov());

  BOOST_CHECK_CLOSE(varZR[0], sigmaLoc1 * sigmaLoc1, 1e-6);
  BOOST_CHECK_SMALL(varZR[1], 1e-12);
}

/// An endcap module measures r and r-phi, so all of the variance lands in r.
/// This pins the scale of the r Jacobian: a factor f in dr/d{x,y} shows as f²
BOOST_AUTO_TEST_CASE(EndcapModule) {
  const double phi = -1.1;
  const double r = 250_mm;
  const Vector3 position(r * std::cos(phi), r * std::sin(phi), 1500_mm);
  const Vector3 radial(std::cos(phi), std::sin(phi), 0);
  const Vector3 rPhi(-std::sin(phi), std::cos(phi), 0);

  auto surface = makePlane(position, radial, rPhi);
  const Vector2 varZR = PixelSpacePointBuilder::computeVarianceZR(
      gctx, *surface, position, localCov());

  BOOST_CHECK_SMALL(varZR[0], 1e-12);
  BOOST_CHECK_CLOSE(varZR[1], sigmaLoc0 * sigmaLoc0, 1e-6);
}

/// A module inclined against both, where the local directions mix into r
BOOST_AUTO_TEST_CASE(InclinedModuleAgainstNumericJacobian) {
  const double phi = 0.7;
  const double r = 180_mm;
  const Vector3 position(r * std::cos(phi), r * std::sin(phi), -800_mm);
  const Vector3 radial(std::cos(phi), std::sin(phi), 0);
  const Vector3 rPhi(-std::sin(phi), std::cos(phi), 0);

  // tilt out of the transverse plane, then rotate within the module plane
  const double tilt = 30 * std::numbers::pi / 180;
  const double skew = 20 * std::numbers::pi / 180;
  const Vector3 tilted =
      std::cos(tilt) * radial + std::sin(tilt) * Vector3::UnitZ();
  const Vector3 loc0 = std::cos(skew) * tilted + std::sin(skew) * rPhi;
  const Vector3 loc1 = -std::sin(skew) * tilted + std::cos(skew) * rPhi;

  auto surface = makePlane(position, loc0, loc1);
  const SquareMatrix2 cov = localCov();
  const Vector2 varZR =
      PixelSpacePointBuilder::computeVarianceZR(gctx, *surface, position, cov);
  const Vector2 expected = numericVarianceZR(*surface, position, cov);

  BOOST_CHECK_CLOSE(varZR[0], expected[0], 1e-3);
  BOOST_CHECK_CLOSE(varZR[1], expected[1], 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
