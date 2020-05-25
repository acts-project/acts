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
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

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

  ///
  /// Coordinate transformations
  ///
  auto jacDirToAngle = detail::jacobianDirectionsToAngles(direction);
  BOOST_TEST(jacDirToAngle != decltype(jacDirToAngle)::Zero());

  auto jacAngleToDir = detail::jacobianAnglesToDirections(direction);
  BOOST_TEST(jacAngleToDir != decltype(jacAngleToDir)::Zero());

  auto product = jacDirToAngle * jacAngleToDir;
  CHECK_CLOSE_ABS(product, decltype(product)::Identity(),
                  std::numeric_limits<double>::epsilon());

  ///
  /// Covariance transport tests
  ///
  // Covariance transport to curvilinear coordinates
  // (1) Local to local
  detail::covarianceTransport(covariance, jacobian, transportJacobian,
                              derivatives, jacobianLocalToGlobal, jacDirToAngle,
                              jacAngleToDir, direction, true);

  // Tests to see that the right components are (un-)changed
  BOOST_REQUIRE_NO_THROW(std::get<BoundSymMatrix>(covariance));
  BOOST_CHECK_NE(std::get<BoundSymMatrix>(covariance) , BoundSymMatrix::Identity());
  BOOST_REQUIRE_NO_THROW(std::get<BoundMatrix>(jacobian));
  BOOST_CHECK_NE(std::get<BoundMatrix>(jacobian) , BoundMatrix(2. * BoundMatrix::Identity()));
  BOOST_CHECK_EQUAL(transportJacobian, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(derivatives, FreeVector::Zero());
  BOOST_CHECK_NE(*jacobianLocalToGlobal, 4. * BoundToFreeMatrix::Identity());
  BOOST_CHECK_EQUAL(direction,
             Vector3D(sqrt(5. / 22.), 3. * sqrt(2. / 55.), 7. / sqrt(110.)));

  // Reset
  covariance = FreeSymMatrix(FreeSymMatrix::Identity());
  jacobian = FreeMatrix(2. * FreeMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  jacobianLocalToGlobal = std::nullopt;

  // (2) Free to local
  detail::covarianceTransport(covariance, jacobian, transportJacobian,
                              derivatives, jacobianLocalToGlobal, jacDirToAngle,
                              jacAngleToDir, direction, true);

  // Tests to see that the right components are (un-)changed
  BOOST_REQUIRE_NO_THROW(std::get<BoundSymMatrix>(covariance));
  BOOST_TEST(std::get<BoundSymMatrix>(covariance) !=
             BoundSymMatrix::Identity());
  BOOST_REQUIRE_NO_THROW(std::get<FreeToBoundMatrix>(jacobian));
  BOOST_TEST(std::get<FreeToBoundMatrix>(jacobian) !=
             FreeToBoundMatrix(2. * FreeToBoundMatrix::Identity()));
  BOOST_TEST(transportJacobian == FreeMatrix::Identity());
  BOOST_TEST(derivatives == FreeVector::Zero());
  BOOST_TEST(jacobianLocalToGlobal.has_value());
  BOOST_TEST(direction ==
             Vector3D(sqrt(5. / 22.), 3. * sqrt(2. / 55.), 7. / sqrt(110.)));

  // Reset
  covariance = FreeSymMatrix(FreeSymMatrix::Identity());
  jacobian = FreeMatrix(2. * FreeMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  jacobianLocalToGlobal = std::nullopt;

  // Repeat transport to free
  // (1) Free to free
  detail::covarianceTransport(covariance, jacobian, transportJacobian,
                              derivatives, jacobianLocalToGlobal, jacDirToAngle,
                              jacAngleToDir, direction, false);

  // Tests to see that the right components are (un-)changed
  BOOST_REQUIRE_NO_THROW(std::get<FreeSymMatrix>(covariance));
  BOOST_TEST(std::get<FreeSymMatrix>(covariance) != FreeSymMatrix::Identity());
  BOOST_REQUIRE_NO_THROW(std::get<FreeMatrix>(jacobian));
  BOOST_TEST(std::get<FreeMatrix>(jacobian) !=
             FreeMatrix(2. * FreeMatrix::Identity()));
  BOOST_TEST(transportJacobian == FreeMatrix::Identity());
  BOOST_TEST(derivatives == FreeVector::Zero());
  BOOST_TEST(!jacobianLocalToGlobal.has_value());
  BOOST_TEST(direction ==
             Vector3D(sqrt(5. / 22.), 3. * sqrt(2. / 55.), 7. / sqrt(110.)));

  // Reset
  covariance = BoundSymMatrix(BoundSymMatrix::Identity());
  jacobian = BoundMatrix(2. * BoundMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  jacobianLocalToGlobal = 4. * BoundToFreeMatrix::Identity();

  // (2) Local to free
  detail::covarianceTransport(covariance, jacobian, transportJacobian,
                              derivatives, jacobianLocalToGlobal, jacDirToAngle,
                              jacAngleToDir, direction, false);

  // Tests to see that the right components are (un-)changed
  BOOST_REQUIRE_NO_THROW(std::get<FreeSymMatrix>(covariance));
  BOOST_TEST(std::get<FreeSymMatrix>(covariance) != FreeSymMatrix::Identity());
  BOOST_REQUIRE_NO_THROW(std::get<BoundToFreeMatrix>(jacobian));
  BOOST_TEST(std::get<BoundToFreeMatrix>(jacobian) !=
             BoundToFreeMatrix(2. * BoundToFreeMatrix::Identity()));
  BOOST_TEST(transportJacobian == FreeMatrix::Identity());
  BOOST_TEST(derivatives == FreeVector::Zero());
  BOOST_TEST(!jacobianLocalToGlobal.has_value());
  BOOST_TEST(direction ==
             Vector3D(sqrt(5. / 22.), 3. * sqrt(2. / 55.), 7. / sqrt(110.)));

  // Reset
  covariance = BoundSymMatrix(BoundSymMatrix::Identity());
  jacobian = BoundMatrix(2. * BoundMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  jacobianLocalToGlobal = 4. * BoundToFreeMatrix::Identity();

  // Repeat transport to surface
  // (1) Local to local
  auto surface = Surface::makeShared<PlaneSurface>(position, direction);
  detail::covarianceTransport(tgContext, covariance, jacobian,
                              transportJacobian, derivatives,
                              jacobianLocalToGlobal, jacDirToAngle,
                              jacAngleToDir, parameters, *surface);
  BOOST_REQUIRE_NO_THROW(std::get<BoundSymMatrix>(covariance));
  BOOST_CHECK_NE(std::get<BoundSymMatrix>(covariance) , BoundSymMatrix::Identity());
  BOOST_REQUIRE_NO_THROW(std::get<BoundMatrix>(jacobian));
  BOOST_CHECK_NE(std::get<BoundMatrix>(jacobian) ,BoundMatrix(2. * BoundMatrix::Identity()));
  BOOST_CHECK_EQUAL(transportJacobian , FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(derivatives, FreeVector::Zero());
  BOOST_CHECK_NE(*jacobianLocalToGlobal ,4. * BoundToFreeMatrix::Identity());
  BOOST_CHECK_EQUAL(parameters , startParameters);

  // (2) Free to local
  covariance = FreeSymMatrix(FreeSymMatrix::Identity());
  jacobian = FreeMatrix(2. * FreeMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  jacobianLocalToGlobal = std::nullopt;
  detail::covarianceTransport(tgContext, covariance, jacobian,
                              transportJacobian, derivatives,
                              jacobianLocalToGlobal, jacDirToAngle,
                              jacAngleToDir, parameters, *surface);

  BOOST_REQUIRE_NO_THROW(std::get<BoundSymMatrix>(covariance));
  BOOST_TEST(std::get<BoundSymMatrix>(covariance) !=
             BoundSymMatrix::Identity());
  BOOST_REQUIRE_NO_THROW(std::get<FreeToBoundMatrix>(jacobian));
  BOOST_TEST(std::get<FreeToBoundMatrix>(jacobian) !=
             FreeToBoundMatrix(2. * FreeToBoundMatrix::Identity()));
  BOOST_TEST(transportJacobian == FreeMatrix::Identity());
  BOOST_TEST(derivatives == FreeVector::Zero());
  BOOST_TEST(jacobianLocalToGlobal.has_value());
  BOOST_TEST(parameters == startParameters);

  ///
  /// State construction tests
  ///
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
  jacobianLocalToGlobal = 4. * BoundToFreeMatrix::Identity();

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
  jacobianLocalToGlobal = 4. * BoundToFreeMatrix::Identity();

  // Produce a bound state with covariance matrix
  boundResult =
      detail::boundState(tgContext, covariance, jacobian, transportJacobian,
                         derivatives, jacobianLocalToGlobal, jacDirToAngle,
                         jacAngleToDir, parameters, true, 1337., *surface);
  BOOST_CHECK(std::get<0>(boundResult).covariance().has_value());
  BOOST_CHECK_NE(*(std::get<0>(boundResult).covariance()),
             BoundSymMatrix(BoundSymMatrix::Identity()));
  BOOST_REQUIRE_NO_THROW(std::get<BoundMatrix>(std::get<1>(boundResult)));
  BOOST_CHECK_NE(std::get<BoundMatrix>(std::get<1>(boundResult)) ,BoundMatrix(2. * BoundMatrix::Identity()));
  BOOST_CHECK_EQUAL(std::get<2>(boundResult), 1337.);
  
  // Produce a free state without covariance matrix
  auto freeResult =
      detail::freeState(covariance, jacobian, transportJacobian,
                         derivatives, jacobianLocalToGlobal, jacDirToAngle,
                         jacAngleToDir, parameters, false, 1337.);
  BOOST_CHECK(!std::get<0>(freeResult).covariance().has_value());
  BOOST_CHECK_EQUAL(std::get<2>(freeResult) , 1337.);

  // Reset
  covariance = FreeSymMatrix(FreeSymMatrix::Identity());
  jacobian = FreeMatrix(2. * FreeMatrix::Identity());
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  jacobianLocalToGlobal = std::nullopt;

  // Produce a bound state with covariance matrix
  freeResult =
      detail::freeState(covariance, jacobian, transportJacobian,
                         derivatives, jacobianLocalToGlobal, jacDirToAngle,
                         jacAngleToDir, parameters, true, 1337.);
  BOOST_CHECK(std::get<0>(freeResult).covariance().has_value());
  BOOST_CHECK_NE(*(std::get<0>(freeResult).covariance()) ,
             FreeSymMatrix(FreeSymMatrix::Identity()));
  BOOST_REQUIRE_NO_THROW(std::get<FreeMatrix>(std::get<1>(freeResult)));
  BOOST_CHECK_NE(std::get<FreeMatrix>(std::get<1>(freeResult)),
             FreeMatrix(2. * FreeMatrix::Identity()));
  BOOST_CHECK_EQUAL(std::get<2>(freeResult) , 1337.);
}
}  // namespace Test
}  // namespace Acts