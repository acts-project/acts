// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/JacobianHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <limits>
#include <numbers>

using namespace Acts;

namespace {
constexpr double eps = 1e-6;

// d(p_x, p_y, p_z) / d(phi, theta, qOverP) evaluated by central finite
// differences of the explicit spherical-to-Cartesian momentum map.
Matrix<3, 3> numericSphericalToFreeMomentumJacobian(double phi, double theta,
                                                     double qOverP,
                                                     double charge) {
  auto momentum = [&](double phi_, double theta_, double qOverP_) {
    return (charge / qOverP_) * makeDirectionFromPhiTheta(phi_, theta_);
  };

  const double h = 1e-6;
  Matrix<3, 3> jacobian;
  jacobian.col(0) =
      (momentum(phi + h, theta, qOverP) - momentum(phi - h, theta, qOverP)) /
      (2 * h);
  jacobian.col(1) =
      (momentum(phi, theta + h, qOverP) - momentum(phi, theta - h, qOverP)) /
      (2 * h);
  jacobian.col(2) =
      (momentum(phi, theta, qOverP + h) - momentum(phi, theta, qOverP - h)) /
      (2 * h);
  return jacobian;
}
}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JacobianHelpersSuite)

BOOST_AUTO_TEST_CASE(SphericalToFreeMomentumJacobianMatchesFiniteDifference) {
  const double phi = 0.6;
  const double theta = 1.1;
  const double charge = 1.;
  const double qOverP = 0.5;
  const double momentum = std::abs(charge / qOverP);

  const Vector3 direction = makeDirectionFromPhiTheta(phi, theta);

  const Matrix<3, 3> analytic =
      sphericalToFreeMomentumJacobian(direction, qOverP, momentum);
  const Matrix<3, 3> numeric =
      numericSphericalToFreeMomentumJacobian(phi, theta, qOverP, charge);

  CHECK_CLOSE_ABS(analytic, numeric, eps);
}

BOOST_AUTO_TEST_CASE(FreeToSphericalMomentumJacobianIsInverse) {
  const double phi = -1.2;
  const double theta = 0.8;
  const double charge = -1.;
  const double qOverP = -0.25;
  const double momentum = std::abs(charge / qOverP);

  const Vector3 direction = makeDirectionFromPhiTheta(phi, theta);
  const Vector3 momentumVector = momentum * direction;

  const Matrix<3, 3> forward =
      sphericalToFreeMomentumJacobian(direction, qOverP, momentum);
  const Matrix<3, 3> inverse =
      freeToSphericalMomentumJacobian(momentumVector, charge);

  CHECK_CLOSE_ABS(inverse * forward, Matrix<3, 3>::Identity(), eps);
  CHECK_CLOSE_ABS(forward * inverse, Matrix<3, 3>::Identity(), eps);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
