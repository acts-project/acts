// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>

namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

using Covariance = BoundSquareMatrix;
using Jacobian = BoundMatrix;

/// These tests do not test for a correct covariance transport but only for the
/// correct conservation or modification of certain variables. A test suite for
/// the numerical correctness is performed in the integration tests.
BOOST_AUTO_TEST_CASE(covariance_engine_test) {
  // Create a test context
  GeometryContext tgContext = GeometryContext();

  auto particleHypothesis = ParticleHypothesis::pion();

  // Build a start vector
  Vector3 position{1., 2., 3.};
  double time = 4.;
  Vector3 direction{sqrt(5. / 22.), 3. * sqrt(2. / 55.), 7. / sqrt(110.)};
  double qop = 0.125;
  FreeVector parameters, startParameters;
  parameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;
  startParameters = parameters;

  // Build covariance matrix, jacobians and related components
  Covariance covariance = Covariance::Identity();
  Jacobian jacobian = 2. * Jacobian::Identity();
  FreeMatrix transportJacobian = 3. * FreeMatrix::Identity();
  FreeVector derivatives;
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  BoundToFreeMatrix boundToFreeJacobian = 4. * BoundToFreeMatrix::Identity();

  // Covariance transport to curvilinear coordinates
  detail::transportCovarianceToCurvilinear(covariance, jacobian,
                                           transportJacobian, derivatives,
                                           boundToFreeJacobian, direction);

  // Tests to see that the right components are (un-)changed
  BOOST_CHECK_NE(covariance, Covariance::Identity());
  BOOST_CHECK_NE(jacobian, 2. * Jacobian::Identity());
  BOOST_CHECK_EQUAL(transportJacobian, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(derivatives, FreeVector::Zero());
  BOOST_CHECK_NE(boundToFreeJacobian, 4. * BoundToFreeMatrix::Identity());
  BOOST_CHECK_EQUAL(
      direction, Vector3(sqrt(5. / 22.), 3. * sqrt(2. / 55.), 7. / sqrt(110.)));

  // Reset
  covariance = Covariance::Identity();
  jacobian = 2. * Jacobian::Identity();
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  boundToFreeJacobian = 4. * BoundToFreeMatrix::Identity();

  // Repeat transport to surface
  FreeToBoundCorrection freeToBoundCorrection(false);
  auto surface = Surface::makeShared<PlaneSurface>(position, direction);
  detail::transportCovarianceToBound(
      tgContext, covariance, jacobian, transportJacobian, derivatives,
      boundToFreeJacobian, parameters, *surface, freeToBoundCorrection);

  BOOST_CHECK_NE(covariance, Covariance::Identity());
  BOOST_CHECK_NE(jacobian, 2. * Jacobian::Identity());
  BOOST_CHECK_EQUAL(transportJacobian, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(derivatives, FreeVector::Zero());
  BOOST_CHECK_NE(boundToFreeJacobian, 4. * BoundToFreeMatrix::Identity());
  BOOST_CHECK_EQUAL(parameters, startParameters);

  // Produce a curvilinear state without covariance matrix
  auto covarianceBefore = covariance;
  auto curvResult = detail::curvilinearState(
      covariance, jacobian, transportJacobian, derivatives, boundToFreeJacobian,
      parameters, particleHypothesis, false, 1337.);
  BOOST_CHECK(std::get<0>(curvResult).covariance().has_value());
  BOOST_CHECK_EQUAL(*(std::get<0>(curvResult).covariance()), covarianceBefore);
  BOOST_CHECK_EQUAL(std::get<2>(curvResult), 1337.);

  // Reset
  covariance = Covariance::Identity();
  jacobian = 2. * Jacobian::Identity();
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  boundToFreeJacobian = 4. * BoundToFreeMatrix::Identity();

  // Produce a curvilinear state with covariance matrix
  curvResult = detail::curvilinearState(
      covariance, jacobian, transportJacobian, derivatives, boundToFreeJacobian,
      parameters, particleHypothesis, true, 1337.);
  BOOST_CHECK(std::get<0>(curvResult).covariance().has_value());
  BOOST_CHECK_NE(*(std::get<0>(curvResult).covariance()),
                 Covariance::Identity());
  BOOST_CHECK_NE(std::get<1>(curvResult), 2. * Jacobian::Identity());
  BOOST_CHECK_EQUAL(std::get<2>(curvResult), 1337.);

  // Produce a bound state without covariance matrix
  covarianceBefore = covariance;
  auto boundResult =
      detail::boundState(tgContext, covariance, jacobian, transportJacobian,
                         derivatives, boundToFreeJacobian, parameters,
                         particleHypothesis, false, 1337., *surface,
                         freeToBoundCorrection)
          .value();
  BOOST_CHECK(std::get<0>(curvResult).covariance().has_value());
  BOOST_CHECK_EQUAL(*(std::get<0>(curvResult).covariance()), covarianceBefore);
  BOOST_CHECK_EQUAL(std::get<2>(boundResult), 1337.);

  // Reset
  covariance = Covariance::Identity();
  jacobian = 2. * Jacobian::Identity();
  transportJacobian = 3. * FreeMatrix::Identity();
  derivatives << 9., 10., 11., 12., 13., 14., 15., 16.;
  boundToFreeJacobian = 4. * BoundToFreeMatrix::Identity();

  // Produce a bound state with covariance matrix
  boundResult =
      detail::boundState(tgContext, covariance, jacobian, transportJacobian,
                         derivatives, boundToFreeJacobian, parameters,
                         ParticleHypothesis::pion(), true, 1337., *surface,
                         freeToBoundCorrection)
          .value();
  BOOST_CHECK(std::get<0>(boundResult).covariance().has_value());
  BOOST_CHECK_NE(*(std::get<0>(boundResult).covariance()),
                 Covariance::Identity());
  BOOST_CHECK_NE(std::get<1>(boundResult), 2. * Jacobian::Identity());
  BOOST_CHECK_EQUAL(std::get<2>(boundResult), 1337.);

  // Reset
  freeToBoundCorrection.apply = true;

  // Produce a bound state with free to bound correction
  boundResult =
      detail::boundState(tgContext, covariance, jacobian, transportJacobian,
                         derivatives, boundToFreeJacobian, parameters,
                         ParticleHypothesis::pion(), true, 1337., *surface,
                         freeToBoundCorrection)
          .value();
  BOOST_CHECK(std::get<0>(boundResult).covariance().has_value());
  BOOST_CHECK_NE(*(std::get<0>(boundResult).covariance()),
                 Covariance::Identity());
}
}  // namespace Test
}  // namespace Acts
