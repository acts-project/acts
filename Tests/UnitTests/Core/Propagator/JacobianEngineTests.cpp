// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Propagator/detail/SympyJacobianEngine.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts;

namespace ActsTest {

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

BOOST_AUTO_TEST_CASE(CanonicalTransportFiniteDifferences) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  const auto source = Surface::makeShared<PlaneSurface>(Transform3::Identity());
  const BoundVector bound(0.4, -0.7, 0.3, 1.1, 0.02, 2.);
  const FreeVector start = transformBoundToFreeParameters(*source, gctx, bound);
  const Vector3 direction = start.segment<3>(eFreeDir0);
  const auto startJacobian =
      source->boundToFreeJacobian(gctx, start.segment<3>(eFreePos0), direction);

  for (double distance : {-10., 10.}) {
    Transform3 placement = Transform3::Identity();
    placement.translate(start.segment<3>(eFreePos0) + distance * direction);
    placement.rotate(AngleAxis3(0.2, Vector3::UnitY()));
    const auto target = Surface::makeShared<PlaneSurface>(placement);
    const Vector3 normal = placement.linear().col(2);
    const auto propagate = [&](const BoundVector& input) {
      FreeVector free = transformBoundToFreeParameters(*source, gctx, input);
      const Vector3 pos = free.segment<3>(eFreePos0);
      const Vector3 dir = free.segment<3>(eFreeDir0);
      const double path =
          normal.dot(placement.translation() - pos) / normal.dot(dir);
      free.segment<3>(eFreePos0) += path * dir;
      return free;
    };
    FreeVector end = propagate(bound);
    FreeMatrix transport = FreeMatrix::Identity();
    transport.block<3, 3>(eFreePos0, eFreeDir0) =
        distance * SquareMatrix3::Identity();
    FreeVector pathDerivative = FreeVector::Zero();
    pathDerivative.segment<3>(eFreePos0) = direction;
    BoundMatrix analytic;
    FreeToBoundMatrix projection;
    BOOST_REQUIRE(detail::boundToBoundTransportJacobian(
                      gctx, *target, end, startJacobian, transport, projection,
                      pathDerivative, analytic)
                      .ok());

    BoundMatrix numerical;
    constexpr double h = 1e-6;
    for (int i = 0; i < static_cast<int>(eBoundSize); ++i) {
      BoundVector plus = bound;
      BoundVector minus = bound;
      plus[i] += h;
      minus[i] -= h;
      const auto upper =
          transformFreeToBoundParameters(propagate(plus), *target, gctx);
      const auto lower =
          transformFreeToBoundParameters(propagate(minus), *target, gctx);
      BOOST_REQUIRE(upper.ok());
      BOOST_REQUIRE(lower.ok());
      numerical.col(i) =
          target->boundParametrization(gctx).difference(*upper, *lower) /
          (2. * h);
    }
    CHECK_CLOSE_ABS(analytic, numerical, 1e-8);

    BoundMatrix sympyJacobian;
    FreeToBoundMatrix sympyProjection;
    BOOST_REQUIRE(detail::sympy::boundToBoundTransportJacobian(
                      gctx, *target, end, transport * startJacobian,
                      sympyProjection, pathDerivative, sympyJacobian)
                      .ok());
    CHECK_CLOSE_ABS(sympyJacobian, numerical, 1e-8);

    BoundMatrix factor = BoundMatrix::Identity();
    factor(eBoundLoc0, eBoundDir0) = 0.2;
    factor(eBoundTime, eBoundDir1) = -0.1;
    const BoundMatrix originalCovariance = factor * factor.transpose();
    BoundMatrix covariance = originalCovariance;
    auto anchorJacobian = startJacobian;
    BOOST_REQUIRE(detail::transportCovarianceToBound(
                      gctx, *target, covariance, analytic, transport,
                      pathDerivative, anchorJacobian, std::nullopt, end,
                      FreeToBoundCorrection{})
                      .ok());
    CHECK_CLOSE_ABS(
        covariance,
        (numerical * originalCovariance * numerical.transpose()).eval(), 1e-6);
    CHECK_CLOSE_ABS(
        anchorJacobian,
        target->boundToFreeJacobian(gctx, end.segment<3>(eFreePos0), direction),
        1e-14);
    BOOST_CHECK_EQUAL(transport, FreeMatrix::Identity());
    BOOST_CHECK_EQUAL(pathDerivative, FreeVector::Zero());
  }
}

/// These tests do not test for a correct covariance transport but only for the
/// correct conservation or modification of certain variables. A test suite for
/// the numerical correctness is performed in the integration tests.
///
BOOST_AUTO_TEST_CASE(jacobian_engine_to_bound) {
  // Create a test context
  GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();

  // Build a start vector
  Vector3 position{1., 2., 3.};
  double time = 4.;
  Vector3 direction = Vector3(5., 2., 7.).normalized();
  double qop = 0.125;

  // Build a surface
  std::shared_ptr<PlaneSurface> pSurface =
      CurvilinearSurface(position, direction).planeSurface();

  // Other rotated surface
  Vector3 odirection = Vector3(6., 2., 8.).normalized();
  std::shared_ptr<PlaneSurface> oSurface =
      CurvilinearSurface(position, odirection).planeSurface();

  // The free parameter vector
  FreeVector freeParameters;
  freeParameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;
  // And its associated bound vector
  BoundVector boundParameters;
  boundParameters << 0., 0., VectorHelpers::phi(direction),
      VectorHelpers::theta(direction), qop, time;

  // Build covariance matrices for bound and free case
  BoundMatrix boundCovariance = 0.025 * BoundMatrix::Identity();
  FreeMatrix freeCovariance = 0.025 * FreeMatrix::Identity();

  FreeMatrix noTransportJacobian = FreeMatrix::Identity();
  FreeMatrix realTransportJacobian = 2 * FreeMatrix::Identity();

  FreeToPathMatrix freeToPathDerivatives =
      pSurface->freeToPathDerivative(tgContext, position, direction);
  BoundToFreeMatrix boundToFreeJacobian =
      pSurface->boundToFreeJacobian(tgContext, position, direction);

  // (1) curvilinear/bound to bound transport jacobian
  // a) test without actual transport into the same surface
  BoundMatrix b2bTransportJacobian;
  FreeToBoundMatrix freeToBoundJacobian;
  BOOST_REQUIRE(detail::boundToBoundTransportJacobian(
                    tgContext, *pSurface, freeParameters, boundToFreeJacobian,
                    noTransportJacobian, freeToBoundJacobian,
                    freeToPathDerivatives, b2bTransportJacobian)
                    .ok());
  BoundMatrix newBoundCovariance =
      b2bTransportJacobian * boundCovariance * b2bTransportJacobian.transpose();
  BOOST_CHECK(boundCovariance.isApprox(newBoundCovariance));
  // b) test without actual transport but to a new surface
  BOOST_REQUIRE(detail::boundToBoundTransportJacobian(
                    tgContext, *oSurface, freeParameters, boundToFreeJacobian,
                    noTransportJacobian, freeToBoundJacobian,
                    freeToPathDerivatives, b2bTransportJacobian)
                    .ok());
  newBoundCovariance =
      b2bTransportJacobian * boundCovariance * b2bTransportJacobian.transpose();
  BOOST_CHECK(!boundCovariance.isApprox(newBoundCovariance));
  // c) test to "the same" surface with transport
  // (not really senseful, but should give a different result)
  BOOST_REQUIRE(detail::boundToBoundTransportJacobian(
                    tgContext, *pSurface, freeParameters, boundToFreeJacobian,
                    realTransportJacobian, freeToBoundJacobian,
                    freeToPathDerivatives, b2bTransportJacobian)
                    .ok());
  newBoundCovariance =
      b2bTransportJacobian * boundCovariance * b2bTransportJacobian.transpose();
  BOOST_CHECK(!boundCovariance.isApprox(newBoundCovariance));

  FreeToBoundMatrix f2bTransportJacobian;
  BOOST_REQUIRE(detail::freeToBoundTransportJacobian(
                    tgContext, *pSurface, freeParameters, noTransportJacobian,
                    freeToPathDerivatives, f2bTransportJacobian)
                    .ok());

  newBoundCovariance =
      f2bTransportJacobian * freeCovariance * f2bTransportJacobian.transpose();
  BOOST_CHECK(!boundCovariance.isApprox(newBoundCovariance));
}

/// These tests do not test for a correct covariance transport but only for the
/// correct conservation or modification of certain variables. A test suite for
/// the numerical correctness is performed in the integration tests.
///
BOOST_AUTO_TEST_CASE(jacobian_engine_to_curvilinear) {
  // Create a test context
  GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();

  // Build a start vector
  Vector3 position{1., 2., 3.};
  Vector3 direction = Vector3(5., 2., 7.).normalized();

  // Build a surface, starting surface for curvilinear
  std::shared_ptr<PlaneSurface> pSurface =
      CurvilinearSurface(position, direction).planeSurface();

  // Build covariance matrices for bound and free case
  BoundMatrix boundCovariance = 0.025 * BoundMatrix::Identity();
  FreeMatrix freeCovariance = 0.025 * FreeMatrix::Identity();

  FreeMatrix noTransportJacobian = FreeMatrix::Identity();
  FreeToBoundMatrix freeToBoundJacobian;

  FreeToPathMatrix freeToPathDerivatives =
      pSurface->freeToPathDerivative(tgContext, position, direction);
  BoundToFreeMatrix boundToFreeJacobian =
      CurvilinearSurface(direction).boundToFreeJacobian();

  // (1) curvilinear/bound to curvilinear transport jacobian
  // a) test without actual transport into the same surface
  BoundMatrix b2cTransportJacobian;
  BOOST_REQUIRE(detail::boundToCurvilinearTransportJacobian(
                    direction, boundToFreeJacobian, noTransportJacobian,
                    freeToBoundJacobian, freeToPathDerivatives,
                    b2cTransportJacobian)
                    .ok());
  BoundMatrix newBoundCovariance =
      b2cTransportJacobian * boundCovariance * b2cTransportJacobian.transpose();
  BOOST_CHECK(boundCovariance.isApprox(newBoundCovariance));
  // b) test to another curvilinear frame at the same point (no transport)
  BOOST_REQUIRE(detail::boundToCurvilinearTransportJacobian(
                    Vector3(4., 5., 6.).normalized(), boundToFreeJacobian,
                    noTransportJacobian, freeToBoundJacobian,
                    freeToPathDerivatives, b2cTransportJacobian)
                    .ok());
  newBoundCovariance =
      b2cTransportJacobian * boundCovariance * b2cTransportJacobian.transpose();
  BOOST_CHECK(!boundCovariance.isApprox(newBoundCovariance));

  FreeToBoundMatrix f2cTransportJacobian =
      detail::freeToCurvilinearTransportJacobian(direction, noTransportJacobian,
                                                 freeToPathDerivatives)
          .value();

  newBoundCovariance =
      f2cTransportJacobian * freeCovariance * f2cTransportJacobian.transpose();
  BOOST_CHECK(!boundCovariance.isApprox(newBoundCovariance));
}

/// These tests do not test for a correct covariance transport but only for the
/// correct conservation or modification of certain variables. A test suite for
/// the numerical correctness is performed in the integration tests.
///
BOOST_AUTO_TEST_CASE(jacobian_engine_to_free) {
  // Create a test context
  GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();

  // Build a start vector
  Vector3 position{1., 2., 3.};
  Vector3 direction = Vector3(5., 2., 7.).normalized();

  // Build a surface, starting surface for curvilinear
  std::shared_ptr<PlaneSurface> pSurface =
      CurvilinearSurface(position, direction).planeSurface();

  // Build covariance matrices for bound and free case
  BoundMatrix boundCovariance = 0.025 * BoundMatrix::Identity();
  FreeMatrix freeCovariance = 0.025 * FreeMatrix::Identity();

  FreeMatrix noTransportJacobian = FreeMatrix::Identity();

  BoundToFreeMatrix boundToFreeJacobian =
      pSurface->boundToFreeJacobian(tgContext, position, direction);

  // (1) bound to free
  BoundToFreeMatrix b2fTransportJacobian = detail::boundToFreeTransportJacobian(
      boundToFreeJacobian, noTransportJacobian);

  FreeMatrix newFreeCovariance1 =
      b2fTransportJacobian * boundCovariance * b2fTransportJacobian.transpose();
  BOOST_CHECK(!newFreeCovariance1.isApprox(freeCovariance));

  // (2) curvilinear to free
  boundToFreeJacobian = CurvilinearSurface(direction).boundToFreeJacobian();
  BoundToFreeMatrix c2fTransportJacobian = detail::boundToFreeTransportJacobian(
      boundToFreeJacobian, noTransportJacobian);

  FreeMatrix newFreeCovariance2 =
      c2fTransportJacobian * boundCovariance * c2fTransportJacobian.transpose();
  BOOST_CHECK(!newFreeCovariance2.isApprox(freeCovariance));
  // But those should be similar/equal
  BOOST_CHECK(newFreeCovariance1.isApprox(newFreeCovariance2));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTest
