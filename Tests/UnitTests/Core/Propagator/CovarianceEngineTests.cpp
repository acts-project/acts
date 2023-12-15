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
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/CovarianceTransport.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Propagator/detail/JacobianEngine.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <memory>
#include <optional>
#include <random>
#include <tuple>
#include <utility>

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

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

BOOST_AUTO_TEST_CASE(CovarianceConversionTest) {
  GeometryContext gctx;

  std::size_t n = 100;

  std::mt19937 rng{42};

  std::uniform_real_distribution<double> locDist(-5., 5.);
  std::uniform_real_distribution<double> globDist(-50., 50.);
  std::uniform_real_distribution<double> angleDist(-2 * M_PI, 2 * M_PI);
  std::uniform_real_distribution<double> unitDist(0, 1);

  for (std::size_t b = 0; b < 10; b++) {
    const Vector3 bField =
        Vector3(unitDist(rng), unitDist(rng), unitDist(rng)) * 3_T;

    for (std::size_t j = 0; j < 10; j++) {
      Transform3 transformA = Transform3::Identity();
      transformA = AngleAxis3(angleDist(rng), Vector3::UnitX());
      transformA = AngleAxis3(angleDist(rng), Vector3::UnitY());
      transformA = AngleAxis3(angleDist(rng), Vector3::UnitZ());

      transformA.translation() << globDist(rng), globDist(rng), globDist(rng);

      auto planeSurfaceA = Surface::makeShared<PlaneSurface>(transformA);

      auto planeSurfaceB = Surface::makeShared<PlaneSurface>(transformA);

      auto boundToBound =
          [&gctx, &bField](
              const auto& parIn, const auto& covIn, const auto& srfA,
              const auto& srfB) -> std::pair<BoundVector, BoundMatrix> {
        Acts::BoundTrackParameters boundParamIn{
            srfA.getSharedPtr(), parIn, covIn, ParticleHypothesis::pion()};

        auto converted =
            boundToBoundConversion(gctx, boundParamIn, srfB, bField).value();

        return {converted.parameters(), converted.covariance().value()};
      };

      BoundMatrix covA;
      covA.setZero();
      covA.diagonal() << 1, 2, 3, 4, 5, 6;

      for (std::size_t i = 0; i < n; i++) {
        BoundVector parA;
        parA << locDist(rng) * 1_mm, locDist(rng) * 1_mm, M_PI / 4.,
            M_PI_2 * 0.9, -1 / 1_GeV, 5_ns;

        // identical surface, this should work
        auto [parB, covB] =
            boundToBound(parA, covA, *planeSurfaceA, *planeSurfaceB);

        // these should be the same because the plane surface are the same
        CHECK_CLOSE_ABS(parA, parB, 1e-9);
        CHECK_CLOSE_COVARIANCE(covA, covB, 1e-9);

        // now go back
        auto [parA2, covA2] =
            boundToBound(parB, covB, *planeSurfaceB, *planeSurfaceA);
        CHECK_CLOSE_ABS(parA, parA2, 1e-9);
        CHECK_CLOSE_COVARIANCE(covA, covA2, 1e-9);
      }

      for (std::size_t i = 0; i < n; i++) {
        BoundVector parA;
        parA << locDist(rng) * 1_mm, locDist(rng) * 1_mm, M_PI / 4.,
            M_PI_2 * 0.9, -1 / 1_GeV, 5_ns;

        // make plane surface that is rotated around the normal
        double angle = angleDist(rng);
        Transform3 transform;
        transform = planeSurfaceA->transform(gctx).rotation();
        transform = AngleAxis3(angle, planeSurfaceA->normal(gctx)) * transform;
        transform.translation() = planeSurfaceA->transform(gctx).translation();

        planeSurfaceB = Surface::makeShared<PlaneSurface>(transform);

        // sanity check that the normal didn't change
        CHECK_CLOSE_ABS(planeSurfaceA->normal(gctx),
                        planeSurfaceB->normal(gctx), 1e-9);

        auto [parB, covB] =
            boundToBound(parA, covA, *planeSurfaceA, *planeSurfaceB);
        BoundVector exp = parA;
        // loc0 and loc1 are rotated
        exp.head<2>() = Eigen::Rotation2D<double>(-angle) * parA.head<2>();

        CHECK_CLOSE_ABS(exp, parB, 1e-9);

        // now go back
        auto [parA2, covA2] =
            boundToBound(parB, covB, *planeSurfaceB, *planeSurfaceA);
        CHECK_CLOSE_ABS(parA, parA2, 1e-9);
        CHECK_CLOSE_COVARIANCE(covA, covA2, 1e-9);
      }

      for (std::size_t i = 0; i < n; i++) {
        // make plane that is slightly rotated
        double angle = unitDist(rng) * M_PI_2;

        Transform3 transform;
        transform = planeSurfaceA->transform(gctx).rotation();

        // figure out rotation axis along local x
        Vector3 axis =
            planeSurfaceA->transform(gctx).rotation() * Vector3::UnitY();
        transform = AngleAxis3(angle, axis) * transform;

        transform.translation() = planeSurfaceA->transform(gctx).translation();

        planeSurfaceB = Surface::makeShared<PlaneSurface>(transform);

        BoundVector parA;
        // loc 0 must be zero so we're on the intersection of both surfaces.
        parA << 0, locDist(rng) * 1_mm, M_PI / 4., M_PI_2 * 0.9, -1 / 1_GeV,
            5_ns;

        auto [parB, covB] =
            boundToBound(parA, covA, *planeSurfaceA, *planeSurfaceB);

        // now go back
        auto [parA2, covA2] =
            boundToBound(parB, covB, *planeSurfaceB, *planeSurfaceA);
        CHECK_CLOSE_ABS(parA, parA2, 1e-9);
        CHECK_CLOSE_COVARIANCE(covA, covA2, 1e-7);
      }

      for (std::size_t i = 0; i < n; i++) {
        // make plane that is slightly rotated
        double angle = unitDist(rng) * M_PI_2;

        Transform3 transform;
        transform = planeSurfaceA->transform(gctx).rotation();

        Vector3 axis =
            planeSurfaceA->transform(gctx).rotation() * Vector3::UnitX();
        transform = AngleAxis3(angle, axis) * transform;

        transform.translation() = planeSurfaceA->transform(gctx).translation();

        planeSurfaceB = Surface::makeShared<PlaneSurface>(transform);

        BoundVector parA;
        // loc 1 must be zero so we're on the intersection of both surfaces.
        parA << locDist(rng) * 1_mm, 0, M_PI / 4., M_PI_2 * 0.9, -1 / 1_GeV,
            5_ns;

        auto [parB, covB] =
            boundToBound(parA, covA, *planeSurfaceA, *planeSurfaceB);

        // now go back
        auto [parA2, covA2] =
            boundToBound(parB, covB, *planeSurfaceB, *planeSurfaceA);
        CHECK_CLOSE_ABS(parA, parA2, 1e-9);
        // tolerance is a bit higher here
        CHECK_CLOSE_COVARIANCE(covA, covA2, 1e-6);
      }

      for (std::size_t i = 0; i < n; i++) {
        BoundVector parA;
        parA << locDist(rng) * 1_mm, locDist(rng) * 1_mm, M_PI / 4.,
            M_PI_2 * 0.9, -1 / 1_GeV, 5_ns;

        Vector3 global = planeSurfaceA->localToGlobal(gctx, parA.head<2>());

        Transform3 transform;
        transform.setIdentity();
        transform.rotate(AngleAxis3(angleDist(rng), Vector3::UnitX()));
        transform.rotate(AngleAxis3(angleDist(rng), Vector3::UnitY()));
        transform.rotate(AngleAxis3(angleDist(rng), Vector3::UnitZ()));
        transform.translation() = global;

        auto perigee = Surface::makeShared<PerigeeSurface>(transform);

        auto [parB, covB] = boundToBound(parA, covA, *planeSurfaceA, *perigee);

        // now go back
        auto [parA2, covA2] =
            boundToBound(parB, covB, *perigee, *planeSurfaceA);
        CHECK_CLOSE_ABS(parA, parA2, 1e-9);
        CHECK_CLOSE_COVARIANCE(covA, covA2, 1e-9);
      }
    }
  }
}
}  // namespace Test
}  // namespace Acts
