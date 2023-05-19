// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

/// Helper function tests
BOOST_AUTO_TEST_CASE(jacobian_engine_helper) {
  // (1a) Standard test with curvilinear not glazingly close to z axis
  Vector3 direction = Vector3(7., 8., 9.).normalized();
  FreeToBoundMatrix f2cJacobian = detail::freeToCurvilinearJacobian(direction);

  ActsScalar phi = VectorHelpers::phi(direction);
  ActsScalar theta = VectorHelpers::theta(direction);
  ActsScalar sinPhi = std::sin(phi);
  ActsScalar cosPhi = std::cos(phi);
  ActsScalar sinTheta = std::sin(theta);
  ActsScalar cosTheta = std::cos(theta);

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

  // (1b) Test with curvilinear glazingly close to z axis
  direction = Vector3(1., 1., 999.).normalized();
  f2cJacobian = detail::freeToCurvilinearJacobian(direction);

  phi = VectorHelpers::phi(direction);
  theta = VectorHelpers::theta(direction);
  sinPhi = std::sin(phi);
  cosPhi = std::cos(phi);
  sinTheta = std::sin(theta);
  cosTheta = std::cos(theta);

  const ActsScalar c = std::hypot(direction.y(), direction.z());
  const ActsScalar invC = 1. / c;
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc0, eFreePos1), -direction.z() * invC,
                  1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc0, eFreePos2), direction.y() * invC,
                  1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc1, eFreePos0), c, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc1, eFreePos1),
                  -direction.x() * direction.y() * invC, 1e-5);
  CHECK_CLOSE_REL(f2cJacobian(eBoundLoc1, eFreePos2),
                  -direction.x() * direction.z() * invC, 1e-5);

  // (2a) Standard test with curvilinear not glazingly close to z axis
  direction = Vector3(7., 8., 9.).normalized();
  BoundToFreeMatrix c2fJacobian = detail::curvilinearToFreeJacobian(direction);

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

  // (3) Test angle - directio relational jacobians
  direction = Vector3(2., 3., 4.).normalized();

  ActsMatrix<7, 8> d2aJacobian = detail::directionToAnglesJacobian(direction);
  ActsMatrix<8, 7> a2dJacobian = detail::anglesToDirectionJacobian(direction);

  auto roundTrip = a2dJacobian * d2aJacobian;
  BOOST_CHECK_EQUAL(roundTrip.cols(), 8);
  BOOST_CHECK_EQUAL(roundTrip.rows(), 8);
}

/// These tests do not test for a correct covariance transport but only for the
/// correct conservation or modification of certain variables. A test suite for
/// the numerical correctness is performed in the integration tests.
///
BOOST_AUTO_TEST_CASE(jacobian_engine_to_bound) {
  // Create a test context
  GeometryContext tgContext = GeometryContext();

  // Build a start vector
  Vector3 position{1., 2., 3.};
  double time = 4.;
  Vector3 direction = Vector3(5., 2., 7.).normalized();
  double qop = 0.125;

  // Build a surface
  auto pSurface = Surface::makeShared<PlaneSurface>(position, direction);

  // Other rotated surface
  Vector3 odirection = Vector3(6., 2., 8.).normalized();
  auto oSurface = Surface::makeShared<PlaneSurface>(position, odirection);

  // The free parameter vector
  FreeVector freeParameters;
  freeParameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;
  // And its associated bound vector
  BoundVector boundParameters;
  boundParameters << 0., 0., VectorHelpers::phi(direction),
      VectorHelpers::theta(direction), qop, time;

  // Build covariance matrices for bound and free case
  BoundSymMatrix boundCovariance = 0.025 * BoundSymMatrix::Identity();
  FreeSymMatrix freeCovariance = 0.025 * FreeSymMatrix::Identity();

  FreeMatrix noTransportJacobian = FreeMatrix::Identity();
  FreeMatrix realTransportJacobian = 2 * FreeMatrix::Identity();

  FreeToPathMatrix freeToPathDerivatives =
      pSurface->freeToPathDerivative(tgContext, freeParameters);
  BoundToFreeMatrix boundToFreeJacobian =
      pSurface->boundToFreeJacobian(tgContext, boundParameters);

  // (1) curvilinear/bound to bound transport jacobian
  // a) test without actual transport into the same surface
  BoundMatrix b2bTransportJacobian = detail::boundToBoundTransportJacobian(
      tgContext, freeParameters, boundToFreeJacobian, noTransportJacobian,
      freeToPathDerivatives, *pSurface);
  BoundSymMatrix newBoundCovariance =
      b2bTransportJacobian * boundCovariance * b2bTransportJacobian.transpose();
  BOOST_CHECK(boundCovariance.isApprox(newBoundCovariance));
  // b) test without actual transport but to a new surface
  b2bTransportJacobian = detail::boundToBoundTransportJacobian(
      tgContext, freeParameters, boundToFreeJacobian, noTransportJacobian,
      freeToPathDerivatives, *oSurface);
  newBoundCovariance =
      b2bTransportJacobian * boundCovariance * b2bTransportJacobian.transpose();
  BOOST_CHECK(not boundCovariance.isApprox(newBoundCovariance));
  // c) test to "the same" surface with transport
  // (not really senseful, but should give a different result)
  b2bTransportJacobian = detail::boundToBoundTransportJacobian(
      tgContext, freeParameters, boundToFreeJacobian, realTransportJacobian,
      freeToPathDerivatives, *pSurface);
  newBoundCovariance =
      b2bTransportJacobian * boundCovariance * b2bTransportJacobian.transpose();
  BOOST_CHECK(not boundCovariance.isApprox(newBoundCovariance));

  // (2) free to bound transport jacobian
  const ActsMatrix<7, 8>& directionToAnglesJacobian =
      detail::directionToAnglesJacobian(direction);
  const ActsMatrix<8, 7>& anglesToDirectionJacobian =
      detail::anglesToDirectionJacobian(direction);

  FreeToBoundMatrix f2bTransportJacobian = detail::freeToBoundTransportJacobian(
      tgContext, freeParameters, directionToAnglesJacobian,
      anglesToDirectionJacobian, noTransportJacobian, freeToPathDerivatives,
      *pSurface);

  newBoundCovariance =
      f2bTransportJacobian * freeCovariance * f2bTransportJacobian.transpose();
  BOOST_CHECK(not boundCovariance.isApprox(newBoundCovariance));
}

/// These tests do not test for a correct covariance transport but only for the
/// correct conservation or modification of certain variables. A test suite for
/// the numerical correctness is performed in the integration tests.
///
BOOST_AUTO_TEST_CASE(jacobian_engine_to_curvilinear) {
  // Create a test context
  GeometryContext tgContext = GeometryContext();

  // Build a start vector
  Vector3 position{1., 2., 3.};
  double time = 4.;
  Vector3 direction = Vector3(5., 2., 7.).normalized();
  double qop = 0.125;

  // Build a surface, starting surface for curvilinear
  auto pSurface = Surface::makeShared<PlaneSurface>(position, direction);

  // The free parameter vector
  FreeVector freeParameters;
  freeParameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;
  // And its associated bound vector
  BoundVector boundParameters;
  boundParameters << 0., 0., VectorHelpers::phi(direction),
      VectorHelpers::theta(direction), qop, time;

  // Build covariance matrices for bound and free case
  BoundSymMatrix boundCovariance = 0.025 * BoundSymMatrix::Identity();
  FreeSymMatrix freeCovariance = 0.025 * FreeSymMatrix::Identity();

  FreeMatrix noTransportJacobian = FreeMatrix::Identity();

  FreeToPathMatrix freeToPathDerivatives =
      pSurface->freeToPathDerivative(tgContext, freeParameters);
  BoundToFreeMatrix boundToFreeJacobian =
      detail::curvilinearToFreeJacobian(direction);

  // (1) curvilinear/bound to curvilinear transport jacobian
  // a) test without actual transport into the same surface
  BoundMatrix b2cTransportJacobian =
      detail::boundToCurvilinearTransportJacobian(
          direction, boundToFreeJacobian, noTransportJacobian,
          freeToPathDerivatives);
  BoundSymMatrix newBoundCovariance =
      b2cTransportJacobian * boundCovariance * b2cTransportJacobian.transpose();
  BOOST_CHECK(boundCovariance.isApprox(newBoundCovariance));
  // b) test to another curvilinear frame at the same point (no transport)
  b2cTransportJacobian = detail::boundToCurvilinearTransportJacobian(
      Vector3(4., 5., 6.).normalized(), boundToFreeJacobian,
      noTransportJacobian, freeToPathDerivatives);
  newBoundCovariance =
      b2cTransportJacobian * boundCovariance * b2cTransportJacobian.transpose();
  BOOST_CHECK(not boundCovariance.isApprox(newBoundCovariance));

  // (2) free to bound transport jacobian
  const ActsMatrix<7, 8>& directionToAnglesJacobian =
      detail::directionToAnglesJacobian(direction);
  const ActsMatrix<8, 7>& anglesToDirectionJacobian =
      detail::anglesToDirectionJacobian(direction);

  FreeToBoundMatrix f2cTransportJacobian =
      detail::freeToCurvilinearTransportJacobian(
          direction, directionToAnglesJacobian, anglesToDirectionJacobian,
          noTransportJacobian, freeToPathDerivatives);

  newBoundCovariance =
      f2cTransportJacobian * freeCovariance * f2cTransportJacobian.transpose();
  BOOST_CHECK(not boundCovariance.isApprox(newBoundCovariance));
}

/// These tests do not test for a correct covariance transport but only for the
/// correct conservation or modification of certain variables. A test suite for
/// the numerical correctness is performed in the integration tests.
///
BOOST_AUTO_TEST_CASE(jacobian_engine_to_free) {
  // Create a test context
  GeometryContext tgContext = GeometryContext();

  // Build a start vector
  Vector3 position{1., 2., 3.};
  double time = 4.;
  Vector3 direction = Vector3(5., 2., 7.).normalized();
  double qop = 0.125;

  // Build a surface, starting surface for curvilinear
  auto pSurface = Surface::makeShared<PlaneSurface>(position, direction);

  // The free parameter vector
  FreeVector freeParameters;
  freeParameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;
  // And its associated bound vector
  BoundVector boundParameters;
  boundParameters << 0., 0., VectorHelpers::phi(direction),
      VectorHelpers::theta(direction), qop, time;

  // Build covariance matrices for bound and free case
  BoundSymMatrix boundCovariance = 0.025 * BoundSymMatrix::Identity();
  FreeSymMatrix freeCovariance = 0.025 * FreeSymMatrix::Identity();

  FreeMatrix noTransportJacobian = FreeMatrix::Identity();

  BoundToFreeMatrix boundToFreeJacobian =
      pSurface->boundToFreeJacobian(tgContext, boundParameters);

  // (1) bound to free
  BoundToFreeMatrix b2fTransportJacobian = detail::boundToFreeTransportJacobian(
      boundToFreeJacobian, noTransportJacobian);

  FreeMatrix newFreeCovariance1 =
      b2fTransportJacobian * boundCovariance * b2fTransportJacobian.transpose();
  BOOST_CHECK(not newFreeCovariance1.isApprox(freeCovariance));

  // (2) curvilinear to free
  boundToFreeJacobian = detail::curvilinearToFreeJacobian(direction);
  BoundToFreeMatrix c2fTransportJacobian = detail::boundToFreeTransportJacobian(
      boundToFreeJacobian, noTransportJacobian);

  FreeMatrix newFreeCovariance2 =
      c2fTransportJacobian * boundCovariance * c2fTransportJacobian.transpose();
  BOOST_CHECK(not newFreeCovariance2.isApprox(freeCovariance));
  // But thos should be similar/equal
  BOOST_CHECK(newFreeCovariance1.isApprox(newFreeCovariance2));
}

}  // namespace Test
}  // namespace Acts