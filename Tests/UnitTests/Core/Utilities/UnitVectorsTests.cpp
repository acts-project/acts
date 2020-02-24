// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <limits>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

using Acts::Vector3D;

namespace {
constexpr auto eps = std::numeric_limits<double>::epsilon();
}

BOOST_AUTO_TEST_SUITE(UnitVectors)

BOOST_AUTO_TEST_CASE(DirectionPhiEta) {
  using Acts::makeDirectionUnitFromPhiEta;

  // along positive x
  const auto xPos1 = makeDirectionUnitFromPhiEta(0.0, 0.0);
  CHECK_CLOSE_REL(xPos1.norm(), 1, eps);
  CHECK_CLOSE_REL(xPos1.dot(Vector3D(1, 0, 0)), 1, eps);
  const auto xPos2 = makeDirectionUnitFromPhiEta(2 * M_PI, 0.0);
  CHECK_CLOSE_REL(xPos2.norm(), 1, eps);
  CHECK_CLOSE_REL(xPos2.dot(Vector3D(1, 0, 0)), 1, eps);
  // along negative x
  const auto xNeg1 = makeDirectionUnitFromPhiEta(M_PI, 0.0);
  CHECK_CLOSE_REL(xNeg1.norm(), 1, eps);
  CHECK_CLOSE_REL(xNeg1.dot(Vector3D(-1, 0, 0)), 1, eps);
  const auto xNeg2 = makeDirectionUnitFromPhiEta(-M_PI, 0.0);
  CHECK_CLOSE_REL(xNeg2.norm(), 1, eps);
  CHECK_CLOSE_REL(xNeg2.dot(Vector3D(-1, 0, 0)), 1, eps);
  // along positive y
  const auto yPos1 = makeDirectionUnitFromPhiEta(M_PI_2, 0.0);
  CHECK_CLOSE_REL(yPos1.norm(), 1, eps);
  CHECK_CLOSE_REL(yPos1.dot(Vector3D(0, 1, 0)), 1, eps);
  const auto yPos2 = makeDirectionUnitFromPhiEta(-3 * M_PI_2, 0.0);
  CHECK_CLOSE_REL(yPos2.norm(), 1, eps);
  CHECK_CLOSE_REL(yPos2.dot(Vector3D(0, 1, 0)), 1, eps);
  // along negative y
  const auto yNeg1 = makeDirectionUnitFromPhiEta(-M_PI_2, 0.0);
  CHECK_CLOSE_REL(yNeg1.norm(), 1, eps);
  CHECK_CLOSE_REL(yNeg1.dot(Vector3D(0, -1, 0)), 1, eps);
  const auto yNeg2 = makeDirectionUnitFromPhiEta(3 * M_PI_2, 0.0);
  CHECK_CLOSE_REL(yNeg2.norm(), 1, eps);
  CHECK_CLOSE_REL(yNeg2.dot(Vector3D(0, -1, 0)), 1, eps);

  const auto inf = std::numeric_limits<double>::infinity();
  // along positive z
  const auto zPos1 = makeDirectionUnitFromPhiEta(0.0, inf);
  CHECK_CLOSE_REL(zPos1.norm(), 1, eps);
  CHECK_CLOSE_REL(zPos1.dot(Vector3D(0, 0, 1)), 1, eps);
  const auto zPos2 = makeDirectionUnitFromPhiEta(M_PI_2, inf);
  CHECK_CLOSE_REL(zPos2.norm(), 1, eps);
  CHECK_CLOSE_REL(zPos2.dot(Vector3D(0, 0, 1)), 1, eps);
  // along negative z
  const auto zNeg1 = makeDirectionUnitFromPhiEta(0.0, -inf);
  CHECK_CLOSE_REL(zNeg1.norm(), 1, eps);
  CHECK_CLOSE_REL(zNeg1.dot(Vector3D(0, 0, -1)), 1, eps);
  const auto zNeg2 = makeDirectionUnitFromPhiEta(M_PI_2, -inf);
  CHECK_CLOSE_REL(zNeg2.norm(), 1, eps);
  CHECK_CLOSE_REL(zNeg2.dot(Vector3D(0, 0, -1)), 1, eps);

  // mixed direction
  const auto mixed1 = makeDirectionUnitFromPhiEta(M_PI_4, 1.0);
  CHECK_CLOSE_REL(mixed1.norm(), 1, eps);
  CHECK_CLOSE_REL(
      mixed1.dot(Vector3D(1, 1, M_SQRT2 * std::sinh(1.0)).normalized()), 1,
      eps);
  const auto mixed2 = makeDirectionUnitFromPhiEta(M_PI_4, -1.0);
  CHECK_CLOSE_REL(mixed2.norm(), 1, eps);
  CHECK_CLOSE_REL(
      mixed2.dot(Vector3D(1, 1, M_SQRT2 * std::sinh(-1.0)).normalized()), 1,
      eps);
  const auto mixed3 = makeDirectionUnitFromPhiEta(-M_PI_4, -1.0);
  CHECK_CLOSE_REL(mixed3.norm(), 1, eps);
  CHECK_CLOSE_REL(
      mixed3.dot(Vector3D(1, -1, M_SQRT2 * std::sinh(-1.0)).normalized()), 1,
      eps);
}

BOOST_AUTO_TEST_CASE(DirectionPhiTheta) {
  using Acts::makeDirectionUnitFromPhiTheta;

  // along positive x
  const auto xPos1 = makeDirectionUnitFromPhiTheta(0.0, M_PI_2);
  CHECK_CLOSE_REL(xPos1.norm(), 1, eps);
  CHECK_CLOSE_REL(xPos1.dot(Vector3D(1, 0, 0)), 1, eps);
  const auto xPos2 = makeDirectionUnitFromPhiTheta(2 * M_PI, M_PI_2);
  CHECK_CLOSE_REL(xPos2.norm(), 1, eps);
  CHECK_CLOSE_REL(xPos2.dot(Vector3D(1, 0, 0)), 1, eps);
  // along negative x
  const auto xNeg1 = makeDirectionUnitFromPhiTheta(M_PI, M_PI_2);
  CHECK_CLOSE_REL(xNeg1.norm(), 1, eps);
  CHECK_CLOSE_REL(xNeg1.dot(Vector3D(-1, 0, 0)), 1, eps);
  const auto xNeg2 = makeDirectionUnitFromPhiTheta(-M_PI, M_PI_2);
  CHECK_CLOSE_REL(xNeg2.norm(), 1, eps);
  CHECK_CLOSE_REL(xNeg2.dot(Vector3D(-1, 0, 0)), 1, eps);
  // along positive y
  const auto yPos1 = makeDirectionUnitFromPhiTheta(M_PI_2, M_PI_2);
  CHECK_CLOSE_REL(yPos1.norm(), 1, eps);
  CHECK_CLOSE_REL(yPos1.dot(Vector3D(0, 1, 0)), 1, eps);
  const auto yPos2 = makeDirectionUnitFromPhiTheta(-3 * M_PI_2, M_PI_2);
  CHECK_CLOSE_REL(yPos2.norm(), 1, eps);
  CHECK_CLOSE_REL(yPos2.dot(Vector3D(0, 1, 0)), 1, eps);
  // along negative y
  const auto yNeg1 = makeDirectionUnitFromPhiTheta(-M_PI_2, M_PI_2);
  CHECK_CLOSE_REL(yNeg1.norm(), 1, eps);
  CHECK_CLOSE_REL(yNeg1.dot(Vector3D(0, -1, 0)), 1, eps);
  const auto yNeg2 = makeDirectionUnitFromPhiTheta(3 * M_PI_2, M_PI_2);
  CHECK_CLOSE_REL(yNeg2.norm(), 1, eps);
  CHECK_CLOSE_REL(yNeg2.dot(Vector3D(0, -1, 0)), 1, eps);

  // along positive z
  const auto zPos1 = makeDirectionUnitFromPhiTheta(0.0, 0.0);
  CHECK_CLOSE_REL(zPos1.norm(), 1, eps);
  CHECK_CLOSE_REL(zPos1.dot(Vector3D(0, 0, 1)), 1, eps);
  const auto zPos2 = makeDirectionUnitFromPhiTheta(M_PI_2, 0.0);
  CHECK_CLOSE_REL(zPos2.norm(), 1, eps);
  CHECK_CLOSE_REL(zPos2.dot(Vector3D(0, 0, 1)), 1, eps);
  // along negative z
  const auto zNeg1 = makeDirectionUnitFromPhiTheta(0.0, M_PI);
  CHECK_CLOSE_REL(zNeg1.norm(), 1, eps);
  CHECK_CLOSE_REL(zNeg1.dot(Vector3D(0, 0, -1)), 1, eps);
  const auto zNeg2 = makeDirectionUnitFromPhiTheta(M_PI_2, M_PI);
  CHECK_CLOSE_REL(zNeg2.norm(), 1, eps);
  CHECK_CLOSE_REL(zNeg2.dot(Vector3D(0, 0, -1)), 1, eps);

  // mixed direction
  const auto mixed1 = makeDirectionUnitFromPhiTheta(M_PI_4, M_PI_4);
  CHECK_CLOSE_REL(mixed1.norm(), 1, eps);
  CHECK_CLOSE_REL(mixed1.dot(Vector3D(1, 1, M_SQRT2).normalized()), 1, eps);
  const auto mixed2 = makeDirectionUnitFromPhiTheta(M_PI_4, 3 * M_PI_4);
  CHECK_CLOSE_REL(mixed2.norm(), 1, eps);
  CHECK_CLOSE_REL(mixed2.dot(Vector3D(1, 1, -M_SQRT2).normalized()), 1, eps);
  const auto mixed3 = makeDirectionUnitFromPhiTheta(-M_PI_4, 3 * M_PI_4);
  CHECK_CLOSE_REL(mixed3.norm(), 1, eps);
  CHECK_CLOSE_REL(mixed3.dot(Vector3D(1, -1, -M_SQRT2).normalized()), 1, eps);
}

namespace {
template <typename Direction, typename RefUnitU, typename RefUnitV>
void testCurvilinear(const Eigen::MatrixBase<Direction>& direction,
                     const Eigen::MatrixBase<RefUnitU>& refU,
                     const Eigen::MatrixBase<RefUnitV>& refV) {
  const auto u = Acts::makeCurvilinearUnitU(direction);
  const auto uv = Acts::makeCurvilinearUnitVectors(direction);
  // verify normalization
  CHECK_CLOSE_ABS(u.norm(), 1, eps);
  CHECK_CLOSE_ABS(uv.first.norm(), 1, eps);
  CHECK_CLOSE_ABS(uv.second.norm(), 1, eps);
  // verify orthonormality
  CHECK_SMALL(u.dot(direction), eps);
  CHECK_SMALL(uv.first.dot(direction), eps);
  CHECK_SMALL(uv.second.dot(direction), eps);
  CHECK_SMALL(uv.first.dot(uv.second), eps);
  // verify u is in the x-y plane
  CHECK_SMALL(u[2], eps);
  CHECK_SMALL(uv.first[2], eps);
  // verify references. do not use element-wise comparison to avoid issues with
  // small, epsilon-like differences.
  CHECK_CLOSE_ABS(u.dot(refU), 1, eps);
  CHECK_CLOSE_ABS(uv.first.dot(refU), 1, eps);
  CHECK_CLOSE_ABS(uv.second.dot(refV), 1, eps);
}
}  // namespace

BOOST_AUTO_TEST_CASE(CurvilinearTransverse) {
  // curvilinear system w/ direction in the transverse plane
  testCurvilinear(Vector3D(1, 1, 0), Vector3D(-1, 1, 0).normalized(),
                  Vector3D(0, 0, 1).normalized());
}

BOOST_AUTO_TEST_CASE(CurvilinearPositiveZ) {
  // curvilinear system w/ direction along z
  testCurvilinear(Vector3D(0, 0, 1), Vector3D(1, 0, 0), Vector3D(0, 1, 0));
}

BOOST_AUTO_TEST_CASE(CurvilinearNegativeZ) {
  // curvilinear system w/ direction along z
  testCurvilinear(Vector3D(0, 0, -1), Vector3D(1, 0, 0), Vector3D(0, -1, 0));
}

BOOST_AUTO_TEST_CASE(CurvilinearCloseToZ) {
  // curvilinear system w/ direction close to z
  testCurvilinear(Vector3D(0, 32 * eps, 1 - 32 * eps), Vector3D(-1, 0, 0),
                  Vector3D(0, -1, 32 * eps));
}

BOOST_AUTO_TEST_SUITE_END()
