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
#include "Acts/SpacePointFormation/StripSpacePointBuilder.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <numbers>

using namespace Acts::UnitLiterals;

namespace Acts::Test {

namespace {

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

// two different strip pitches, so that the two directions cannot be confused
constexpr double var1 = 4.7e-4;
constexpr double var2 = 9.0e-4;

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

/// The two local variances of a strip crossing, from the information matrix of
/// the two measurements built and inverted here
Vector2 referenceVariances(const double theta) {
  const Vector2 n1(1, 0);
  const Vector2 n2(std::cos(theta), std::sin(theta));
  const SquareMatrix2 information =
      n1 * n1.transpose() / var1 + n2 * n2.transpose() / var2;
  return information.inverse().diagonal();
}

/// A barrel module, whose only direction reaching `varianceZ` is the one
/// along z: pointing loc0 and then loc1 at it reads the two out one at a time
Vector2 varianceZR(const bool precisionAlongZ, const double theta) {
  const double phi = 0.4;
  const double r = 400_mm;
  const Vector3 position(r * std::cos(phi), r * std::sin(phi), 300_mm);
  const Vector3 rPhi(-std::sin(phi), std::cos(phi), 0);
  const Vector3 z = Vector3::UnitZ();

  auto surface = precisionAlongZ ? makePlane(position, z, rPhi)
                                 : makePlane(position, rPhi, z);
  return StripSpacePointBuilder::computeVarianceZR(gctx, *surface, position,
                                                   var1, var2, theta);
}

// orthogonal strips and the ITk strip stereo angle
constexpr std::array<double, 3> testAngles{std::numbers::pi / 2, 0.4, 40e-3};

}  // namespace

BOOST_AUTO_TEST_SUITE(StripSpacePointBuilderTests)

/// Strip2 adds nothing along the precision direction of strip1, which keeps
/// its own variance whatever the stereo angle
BOOST_AUTO_TEST_CASE(PrecisionDirection) {
  for (const double theta : testAngles) {
    BOOST_TEST_CONTEXT("theta = " << theta) {
      BOOST_CHECK_CLOSE(varianceZR(true, theta)[0],
                        referenceVariances(theta)[0], 1e-6);
      BOOST_CHECK_CLOSE(varianceZR(true, theta)[0], var1, 1e-6);
    }
  }
}

/// Along the strips the crossing is located to 1 / sin(theta) of a pitch, so
/// this direction degrades as the two become parallel
BOOST_AUTO_TEST_CASE(AlongStripDirection) {
  for (const double theta : testAngles) {
    BOOST_TEST_CONTEXT("theta = " << theta) {
      BOOST_CHECK_CLOSE(varianceZR(false, theta)[0],
                        referenceVariances(theta)[1], 1e-6);
    }
  }
}

/// Orthogonal strips measure the two local directions independently
BOOST_AUTO_TEST_CASE(OrthogonalStrips) {
  const double theta = std::numbers::pi / 2;
  BOOST_CHECK_CLOSE(varianceZR(true, theta)[0], var1, 1e-6);
  BOOST_CHECK_CLOSE(varianceZR(false, theta)[0], var2, 1e-6);
}

/// A barrel module carries no radial information either way
BOOST_AUTO_TEST_CASE(NoRadialVariance) {
  BOOST_CHECK_SMALL(varianceZR(true, 40e-3)[1], 1e-12);
  BOOST_CHECK_SMALL(varianceZR(false, 40e-3)[1], 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
