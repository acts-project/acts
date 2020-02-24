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

BOOST_AUTO_TEST_SUITE(UnitVectors)

BOOST_AUTO_TEST_CASE(DirectionTransverse) {
  // curvilinear system w/ direction in the transverse plane
  testCurvilinear(Vector3D(1, 1, 0), Vector3D(-1, 1, 0).normalized(),
                  Vector3D(0, 0, 1).normalized());
}

BOOST_AUTO_TEST_CASE(DirectionPositiveZ) {
  // curvilinear system w/ direction along z
  testCurvilinear(Vector3D(0, 0, 1), Vector3D(1, 0, 0), Vector3D(0, 1, 0));
}

BOOST_AUTO_TEST_CASE(DirectionNegativeZ) {
  // curvilinear system w/ direction along z
  testCurvilinear(Vector3D(0, 0, -1), Vector3D(1, 0, 0), Vector3D(0, -1, 0));
}

BOOST_AUTO_TEST_CASE(DirectionCloseToZ) {
  // curvilinear system w/ direction close to z
  testCurvilinear(Vector3D(0, 32 * eps, 1 - 32 * eps), Vector3D(-1, 0, 0),
                  Vector3D(0, -1, 32 * eps));
}

BOOST_AUTO_TEST_SUITE_END()
