// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

using Covariance = std::variant<BoundSymMatrix, FreeSymMatrix>;
using Jacobian =
    std::variant<BoundMatrix, FreeToBoundMatrix, BoundToFreeMatrix, FreeMatrix>;

/// These tests do not test for a correct covariance transport but only for the
/// correct conservation or modification of certain variables. A test suite for
/// the numerical correctness is performed in the integration tests.
BOOST_AUTO_TEST_CASE(covariance_engine_test) {
  // Create a test context
  GeometryContext tgContext = GeometryContext();

  // Build a start vector
  Vector3D position{1., 2., 3.};
  double time = 4.;
  Vector3D direction{sqrt(5. / 22.), 3. * sqrt(2. / 55.), 7. / sqrt(110.)};
  double qop = 0.125;
  FreeVector parameters, startParameters;
  parameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;
  startParameters = parameters;

  // Build covariance matrix, jacobians and related components
  Covariance covariance = BoundSymMatrix(BoundSymMatrix::Identity());
  Jacobian jacobian = BoundMatrix(2. * BoundMatrix::Identity());
  FreeMatrix transportJacobian = 3. * FreeMatrix::Identity();
  FreeVector derivatives;
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  std::optional<BoundToFreeMatrix> jacobianLocalToGlobal =
      4. * BoundToFreeMatrix::Identity();
  ActsMatrixD<8, 7> jacDirToAngle = ActsMatrixD<8, 7>::Zero();
  ActsMatrixD<7, 8> jacAngleToDir = ActsMatrixD<7, 8>::Zero();

  // Covariance transport to curvilinear coordinates
  detail::covarianceTransport(covariance, jacobian, transportJacobian,
                              derivatives, jacobianLocalToGlobal, jacDirToAngle,
                              jacAngleToDir, direction, true);

  // Tests to see that the right components are (un-)changed
  BOOST_CHECK_NE(std::get<BoundSymMatrix>(covariance) , BoundSymMatrix::Identity());
  BOOST_CHECK_NE(std::get<BoundMatrix>(jacobian) , BoundMatrix(2. * BoundMatrix::Identity()));
  BOOST_CHECK_EQUAL(transportJacobian, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(derivatives, FreeVector::Zero());
  BOOST_CHECK_NE(*jacobianLocalToGlobal, 4. * BoundToFreeMatrix::Identity());
  BOOST_CHECK_EQUAL(direction,
             Vector3D(sqrt(5. / 22.), 3. * sqrt(2. / 55.), 7. / sqrt(110.)));

  // Reset
  covariance = BoundSymMatrix(BoundSymMatrix::Identity());
  jacobian = BoundMatrix(2. * BoundMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  *jacobianLocalToGlobal = 4. * BoundToFreeMatrix::Identity();

  // Repeat transport to surface
  auto surface = Surface::makeShared<PlaneSurface>(position, direction);
  detail::covarianceTransport(tgContext, covariance, jacobian,
                              transportJacobian, derivatives,
                              jacobianLocalToGlobal, jacDirToAngle,
                              jacAngleToDir, parameters, *surface);

  BOOST_CHECK_NE(std::get<BoundSymMatrix>(covariance) , BoundSymMatrix::Identity());
  BOOST_CHECK_NE(std::get<BoundMatrix>(jacobian) ,BoundMatrix(2. * BoundMatrix::Identity()));
  BOOST_CHECK_EQUAL(transportJacobian , FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(derivatives, FreeVector::Zero());
  BOOST_CHECK_NE(*jacobianLocalToGlobal ,4. * BoundToFreeMatrix::Identity());
  BOOST_CHECK_EQUAL(parameters , startParameters);

  // Produce a curvilinear state without covariance matrix
  auto covarianceBefore = covariance;
  auto curvResult = detail::curvilinearState(
      covariance, jacobian, transportJacobian, derivatives,
      jacobianLocalToGlobal, jacDirToAngle, jacAngleToDir, parameters, false,
      1337.);
  BOOST_CHECK(std::get<0>(curvResult).covariance().has_value());
  BOOST_CHECK_EQUAL(*(std::get<0>(curvResult).covariance()), covarianceBefore);
  BOOST_CHECK_EQUAL(std::get<2>(curvResult), 1337.);

  // Reset
  covariance = BoundSymMatrix(BoundSymMatrix::Identity());
  jacobian = BoundMatrix(2. * BoundMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  *jacobianLocalToGlobal = 4. * BoundToFreeMatrix::Identity();

  // Produce a curvilinear state with covariance matrix
  curvResult = detail::curvilinearState(covariance, jacobian, transportJacobian,
                                        derivatives, jacobianLocalToGlobal,
                                        jacDirToAngle, jacAngleToDir,
                                        parameters, true, 1337.);
  BOOST_CHECK(std::get<0>(curvResult).covariance().has_value());
  BOOST_CHECK_NE(*(std::get<0>(curvResult).covariance()) , BoundSymMatrix::Identity());
  BOOST_CHECK_NE(std::get<BoundMatrix>(std::get<1>(curvResult)) ,BoundMatrix(2. * BoundMatrix::Identity()));
  BOOST_CHECK_EQUAL(std::get<2>(curvResult), 1337.);

  // Produce a bound state without covariance matrix
  covarianceBefore = covariance;
  auto boundResult =
      detail::boundState(tgContext, covariance, jacobian, transportJacobian,
                         derivatives, jacobianLocalToGlobal, jacDirToAngle,
                         jacAngleToDir, parameters, false, 1337., *surface);
  BOOST_CHECK(std::get<0>(boundResult).covariance().has_value());
  BOOST_CHECK_EQUAL(*(std::get<0>(boundResult).covariance()), covarianceBefore);
  BOOST_CHECK_EQUAL(std::get<2>(boundResult), 1337.);

  // Reset
  covariance = BoundSymMatrix(BoundSymMatrix::Identity());
  jacobian = BoundMatrix(2. * BoundMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  *jacobianLocalToGlobal = 4. * BoundToFreeMatrix::Identity();

  // Produce a bound state with covariance matrix
  boundResult =
      detail::boundState(tgContext, covariance, jacobian, transportJacobian,
                         derivatives, jacobianLocalToGlobal, jacDirToAngle,
                         jacAngleToDir, parameters, true, 1337., *surface);
  BOOST_CHECK(std::get<0>(boundResult).covariance().has_value());
  BOOST_CHECK_NE(*(std::get<0>(boundResult).covariance()),
             BoundSymMatrix(BoundSymMatrix::Identity()));
  BOOST_CHECK_NE(std::get<BoundMatrix>(std::get<1>(boundResult)) ,BoundMatrix(2. * BoundMatrix::Identity()));
  BOOST_CHECK_EQUAL(std::get<2>(boundResult), 1337.);
}
}  // namespace Test
}  // namespace Acts