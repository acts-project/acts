// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/BoundParametrization.hpp"
#include "Acts/Utilities/JacobianHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <limits>
#include <numbers>
#include <type_traits>

using namespace Acts;

namespace ActsTests {
BOOST_AUTO_TEST_SUITE(BoundParametrizationSuite)

constexpr auto chart = BoundParametrization::globalSpherical();
static_assert(chart.isGlobalSpherical());
static_assert(std::is_trivially_copyable_v<BoundParametrization>);
static_assert(std::is_nothrow_copy_constructible_v<BoundParametrization>);
static_assert(noexcept(BoundParametrization::globalSpherical()));

BOOST_AUTO_TEST_CASE(ValuesMapsAndExpressions) {
  BoundVector values(1., 2., 0.3, 1.2, 0.01, 4.);
  Eigen::Map<BoundVector> map(values.data());
  const Eigen::Map<const BoundVector> constMap(values.data());
  auto copied = chart;
  const Vector3 expected = makeDirectionFromPhiTheta(0.3, 1.2);
  BOOST_CHECK_EQUAL(copied.direction(map), expected);
  BOOST_CHECK_EQUAL(chart.direction(constMap), expected);
  BOOST_CHECK_EQUAL(chart.direction(values + BoundVector::Zero()), expected);
  CHECK_CLOSE_ABS(chart.directionCoordinates(3. * expected), Vector2(0.3, 1.2),
                  1e-14);
  BOOST_CHECK(chart.isValid(constMap));
  BOOST_CHECK(chart.isInDomain(constMap));
  BOOST_CHECK_EQUAL(chart.canonicalize(map), values);
  BOOST_CHECK_EQUAL(chart.difference(map, constMap), BoundVector::Zero());
  BOOST_CHECK_EQUAL(chart.update(map, BoundVector::Zero()), values);
}

BOOST_AUTO_TEST_CASE(Derivatives) {
  constexpr double h = 1e-6;
  for (double phi : {-2.8, -0.2, 1.4, 3.0}) {
    for (double theta : {0.01, 0.7, 1.5, 3.13}) {
      BoundVector values(1., 2., phi, theta, 0.01, 4.);
      const Vector3 dir = chart.direction(values);
      const auto forward = chart.directionJacobian(values);
      const auto inverse = chart.directionCoordinatesJacobian(dir);
      // The legacy helpers reconstruct sin(theta) from 1 - cos(theta)^2,
      // which loses precision near the poles. Finite differences below test
      // the new derivatives independently.
      CHECK_CLOSE_OR_SMALL(forward, sphericalToFreeDirectionJacobian(dir),
                           1e-10, 1e-12);
      CHECK_CLOSE_OR_SMALL(inverse, freeToSphericalDirectionJacobian(dir),
                           1e-10, 1e-12);
      CHECK_CLOSE_ABS((inverse * forward).eval(), SquareMatrix2::Identity(),
                      1e-12);
      for (int i = 0; i < 2; ++i) {
        auto plus = values;
        auto minus = values;
        plus[eBoundDir0 + i] += h;
        minus[eBoundDir0 + i] -= h;
        const Vector3 numerical =
            (chart.direction(plus) - chart.direction(minus)) / (2. * h);
        CHECK_CLOSE_ABS(forward.col(i), numerical, 1e-8);
      }
      const Vector3 scaled = 2.5 * dir;
      const auto scaledInverse = chart.directionCoordinatesJacobian(scaled);
      for (int i = 0; i < 3; ++i) {
        Vector3 plus = scaled, minus = scaled;
        plus[i] += h;
        minus[i] -= h;
        const Vector2 numerical = (chart.directionCoordinates(plus) -
                                   chart.directionCoordinates(minus)) /
                                  (2. * h);
        CHECK_CLOSE_ABS(scaledInverse.col(i), numerical, 1e-6);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(PolesAndInvalidInputs) {
  BoundVector values = BoundVector::Zero();
  BOOST_CHECK(chart.isValid(values));
  BOOST_CHECK(!chart.isInDomain(values));
  BOOST_CHECK_EQUAL(chart.direction(values), Vector3::UnitZ());
  BOOST_CHECK_EQUAL(chart.directionCoordinates(Vector3::UnitZ()),
                    Vector2::Zero());
  BOOST_CHECK_EQUAL(chart.directionCoordinates(-Vector3::UnitZ()),
                    Vector2(0., std::numbers::pi));
  BOOST_CHECK(chart.directionJacobian(values).allFinite());
  BOOST_CHECK_THROW(chart.directionCoordinatesJacobian(Vector3::UnitZ()),
                    std::domain_error);
  BOOST_CHECK_THROW(chart.canonicalizationJacobian(values), std::domain_error);
  values[eBoundDir1] = std::numbers::pi;
  BOOST_CHECK(!chart.isInDomain(values, 0.));
  BOOST_CHECK_THROW(chart.canonicalizationJacobian(values), std::domain_error);
  values[eBoundDir1] = 1e-8;
  BOOST_CHECK(chart.isInDomain(values, 1e-9));
  BOOST_CHECK(!chart.isInDomain(values, 1e-7));
  BOOST_CHECK_THROW(
      chart.directionCoordinatesJacobian(chart.direction(values), 1e-7),
      std::domain_error);
  BOOST_CHECK_THROW(chart.isInDomain(values, -1.), std::invalid_argument);
  BOOST_CHECK_THROW(chart.directionCoordinates(Vector3::Zero()),
                    std::domain_error);
  values[eBoundDir0] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK(!chart.isValid(values));
  BOOST_CHECK(!chart.isInDomain(values));
  BOOST_CHECK_THROW(chart.direction(values), std::domain_error);
  BOOST_CHECK_THROW(chart.canonicalize(values), std::domain_error);
  BOOST_CHECK_THROW(chart.difference(values, BoundVector::Zero()),
                    std::domain_error);
}

BOOST_AUTO_TEST_CASE(PeriodicResidualAndCovarianceBasis) {
  BoundVector lhs(1., 2., std::numbers::pi - 0.01, 1.2, 0.01, 4.);
  BoundVector rhs = lhs;
  rhs[eBoundDir0] = -std::numbers::pi + 0.01;
  const BoundVector delta = chart.difference(lhs, rhs);
  BOOST_CHECK_CLOSE(delta[eBoundDir0], -0.02, 1e-9);
  CHECK_CLOSE_ABS(chart.update(rhs, delta), lhs, 1e-14);

  BoundMatrix covariance = BoundMatrix::Identity();
  covariance(eBoundLoc0, eBoundDir1) = 0.2;
  covariance(eBoundDir1, eBoundLoc0) = 0.2;
  covariance(eBoundDir0, eBoundDir1) = -0.3;
  covariance(eBoundDir1, eBoundDir0) = -0.3;
  const Eigen::Map<const BoundMatrix> covMap(covariance.data());
  for (double theta : {-0.7, 0.7, 3.8, 7.0}) {
    BoundVector raw = lhs;
    raw[eBoundDir1] = theta;
    const Eigen::Map<const BoundVector> map(raw.data());
    const auto [canonical, transformed] = chart.canonicalize(map, covMap);
    CHECK_CLOSE_ABS(chart.direction(canonical), chart.direction(raw), 1e-14);
    const BoundMatrix jacobian = chart.canonicalizationJacobian(map);
    const double sign = (theta == -0.7 || theta == 3.8) ? -1. : 1.;
    BOOST_CHECK_EQUAL(transformed(eBoundLoc0, eBoundDir1), sign * 0.2);
    BOOST_CHECK_EQUAL(transformed(eBoundDir0, eBoundDir1), sign * -0.3);
    BOOST_CHECK_EQUAL(transformed(eBoundDir1, eBoundDir1), 1.);
    constexpr double h = 1e-6;
    for (int i = 0; i < eBoundSize; ++i) {
      BoundVector plus = raw, minus = raw;
      plus[i] += h;
      minus[i] -= h;
      const BoundVector numerical =
          chart.difference(chart.canonicalize(plus),
                           chart.canonicalize(minus)) /
          (2. * h);
      CHECK_CLOSE_ABS(jacobian.col(i), numerical, 1e-8);
    }
    const Matrix<3, 2> before = chart.directionJacobian(raw);
    const Matrix<3, 2> after = chart.directionJacobian(canonical);
    CHECK_CLOSE_ABS((before * covariance.block<2, 2>(eBoundDir0, eBoundDir0) *
                     before.transpose())
                        .eval(),
                    (after * transformed.block<2, 2>(eBoundDir0, eBoundDir0) *
                     after.transpose())
                        .eval(),
                    1e-13);
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
