// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/Interactions.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_CASE(PhysicalScatteringAndProjection) {
  const MaterialSlab slab(makeSilicon(), 1.);
  const auto particle = ParticleHypothesis::pion();
  const auto effects =
      detail::computeMaterialEffects(slab, particle, 0.1f, true, true, true);
  const double theta0 = computeMultipleScatteringTheta0(
      slab, particle.absolutePdg(), particle.mass(), 0.1f,
      particle.absoluteCharge());
  BOOST_CHECK_EQUAL(effects.scatteringVariance, theta0 * theta0);
  BOOST_CHECK_GT(effects.eLoss, 0.);
  BOOST_CHECK_GT(effects.varianceQoverP, 0.);
  const auto parametrization = BoundParametrization::globalSpherical();
  for (const Vector3& direction :
       {Vector3(1., 2., 3.).normalized().eval(), Vector3::UnitX().eval(),
        Vector3(1e-7, 0., 1.).normalized().eval()}) {
    const auto noise = detail::scatteringNoise(parametrization, direction,
                                               effects.scatteringVariance);
    BOOST_REQUIRE(noise.ok());
    const auto jacobian =
        parametrization.directionCoordinatesJacobianFromUnitDirection(direction,
                                                                      0.);
    BOOST_REQUIRE(jacobian.ok());
    const Matrix<3, 3> cartesian =
        effects.scatteringVariance *
        (Matrix<3, 3>::Identity() - direction * direction.transpose());
    const Matrix<2, 2> projected =
        (*jacobian) * cartesian * jacobian->transpose();
    BOOST_CHECK(noise->isApprox(projected, 1e-12));
    BOOST_CHECK_EQUAL((*noise)(0, 1), 0.);
    BOOST_CHECK_EQUAL((*noise)(1, 1), theta0 * theta0);
  }
  const auto noCovariance =
      detail::computeMaterialEffects(slab, particle, 0.1f, true, true, false);
  BOOST_CHECK_EQUAL(noCovariance.scatteringVariance, 0.);
  BOOST_CHECK_EQUAL(noCovariance.varianceQoverP, 0.);
  BOOST_CHECK_EQUAL(noCovariance.eLoss, effects.eLoss);
}

BOOST_AUTO_TEST_CASE(AxialScatteringDiagnostics) {
  const auto parametrization = BoundParametrization::globalSpherical();
  for (double sign : {-1., 1.}) {
    const Vector3 direction(0., 0., sign);
    const auto noise =
        detail::scatteringNoise(parametrization, direction, 0.01);
    BOOST_REQUIRE(!noise.ok());
    BOOST_CHECK(noise.error() == BoundParametrizationError::SingularDirection);
    const Vector2 diagnostic =
        detail::globalAngularScatteringVariances(direction, 0.01);
    BOOST_CHECK(std::isinf(diagnostic[0]));
    BOOST_CHECK_EQUAL(diagnostic[1], 0.01);
    const auto zero = detail::scatteringNoise(parametrization, direction, 0.);
    BOOST_REQUIRE(zero.ok());
    BOOST_CHECK(zero->isZero());
    BOOST_CHECK(
        detail::globalAngularScatteringVariances(direction, 0.).isZero());
  }
}

BOOST_AUTO_TEST_CASE(MaterialFailurePreservesState) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext mctx;
  const StraightLineStepper stepper;
  const StraightLineStepper::Options options(gctx, mctx);
  struct State {
    struct {
      Direction direction = Direction::Forward();
    } options;
    StraightLineStepper::State stepping;
  };
  State state{{}, stepper.makeState(options)};
  const auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity());
  const BoundMatrix covariance = BoundMatrix::Identity();
  stepper.initialize(
      state.stepping,
      BoundTrackParameters(surface, BoundVector(0., 0., 0.3, 1.2, 0.1, 0.),
                           covariance, ParticleHypothesis::pion()));
  stepper.update(state.stepping, Vector3::Zero(), Vector3::UnitZ(), 0.1, 0.);
  const double momentum = stepper.absoluteMomentum(state.stepping);
  const MaterialSlab slab(makeSilicon(), 1.);
  const auto result = detail::performMaterialInteraction(
      state, stepper, slab, surface->boundParametrization(gctx),
      NoiseUpdateMode::addNoise, true, true);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() == BoundParametrizationError::SingularDirection);
  BOOST_CHECK_EQUAL(stepper.absoluteMomentum(state.stepping), momentum);
  BOOST_CHECK_EQUAL(*stepper.covariance(state.stepping), covariance);
  state.stepping.cov.reset();
  BOOST_CHECK(detail::performMaterialInteraction(
                  state, stepper, slab, surface->boundParametrization(gctx),
                  NoiseUpdateMode::addNoise, true, true)
                  .ok());
  BOOST_CHECK_LT(stepper.absoluteMomentum(state.stepping), momentum);

  // A failed noise subtraction must not commit the energy-loss update either.
  stepper.update(state.stepping, Vector3::Zero(), Vector3::UnitX(), 0.1, 0.);
  state.stepping.cov = BoundMatrix::Zero();
  const double beforeRemoval = stepper.absoluteMomentum(state.stepping);
  const auto failedRemoval = detail::performMaterialInteraction(
      state, stepper, slab, surface->boundParametrization(gctx),
      NoiseUpdateMode::removeNoise, true, true);
  BOOST_REQUIRE(!failedRemoval.ok());
  BOOST_CHECK(failedRemoval.error() ==
              MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval);
  BOOST_CHECK_EQUAL(stepper.absoluteMomentum(state.stepping), beforeRemoval);
  BOOST_CHECK(stepper.covariance(state.stepping)->isZero(0.));

  // Addition and removal retain the canonical diagonal update convention.
  const Vector3 direction = Vector3(3., 4., 5.).normalized();
  stepper.update(state.stepping, Vector3::Zero(), direction, 0.1, 0.);
  state.stepping.cov = covariance;
  const auto added = detail::performMaterialInteraction(
      state, stepper, slab, surface->boundParametrization(gctx),
      NoiseUpdateMode::addNoise, true, false);
  BOOST_REQUIRE(added.ok());
  const auto expected = detail::globalAngularScatteringVariances(
      direction, added->scatteringVariance);
  const BoundMatrix after = *stepper.covariance(state.stepping);
  BOOST_CHECK_CLOSE(after(eBoundDir0, eBoundDir0), 1. + expected[0], 1e-10);
  BOOST_CHECK_CLOSE(after(eBoundDir1, eBoundDir1), 1. + expected[1], 1e-10);
  BOOST_CHECK_EQUAL(after(eBoundQOverP, eBoundQOverP), 1.);
  const auto removed = detail::performMaterialInteraction(
      state, stepper, slab, surface->boundParametrization(gctx),
      NoiseUpdateMode::removeNoise, true, false);
  BOOST_REQUIRE(removed.ok());
  BOOST_CHECK(stepper.covariance(state.stepping)->isApprox(covariance, 1e-12));
}

BOOST_AUTO_TEST_CASE(CorrelatedMaterialNoise) {
  BoundMatrix covariance = BoundMatrix::Identity();
  Matrix<2, 2> noise;
  noise << 0.2, -0.1, -0.1, 0.3;
  BOOST_REQUIRE(detail::updateMaterialCovariance(covariance, noise, 0.4,
                                                 NoiseUpdateMode::addNoise)
                    .ok());
  BOOST_CHECK_EQUAL(covariance(eBoundDir0, eBoundDir1), -0.1);
  BOOST_CHECK_EQUAL(covariance(eBoundDir1, eBoundDir0), -0.1);
  BOOST_REQUIRE(detail::updateMaterialCovariance(covariance, noise, 0.4,
                                                 NoiseUpdateMode::removeNoise)
                    .ok());
  BOOST_CHECK(covariance.isApprox(BoundMatrix::Identity(), 1e-14));

  // A nonzero rank-deficient result is valid as well.
  covariance.setIdentity();
  covariance(eBoundDir0, eBoundDir1) = 0.5;
  covariance(eBoundDir1, eBoundDir0) = 0.5;
  BOOST_REQUIRE(detail::updateMaterialCovariance(
                    covariance, 0.5 * Matrix<2, 2>::Identity(), 0.,
                    NoiseUpdateMode::removeNoise)
                    .ok());
  BOOST_CHECK_EQUAL(covariance(eBoundDir0, eBoundDir0), 0.5);

  // Exact removal at a semidefinite boundary is allowed.
  covariance.setZero();
  covariance.block<2, 2>(eBoundDir0, eBoundDir0) = noise;
  BOOST_REQUIRE(detail::updateMaterialCovariance(covariance, noise, 0.,
                                                 NoiseUpdateMode::removeNoise)
                    .ok());
  BOOST_CHECK(covariance.isZero(0.));
}

BOOST_AUTO_TEST_CASE(NoiseRemovalChecksAllCorrelations) {
  Matrix<2, 2> noise = Matrix<2, 2>::Zero();
  noise(0, 0) = 0.5;
  // Positive direction diagonals alone do not establish positivity of the
  // complete covariance: the position/direction correlation also matters.
  for (double unitScale : {1e-150, 1e-12, 1., 1e12, 1e150}) {
    BoundMatrix covariance = BoundMatrix::Identity();
    covariance(eBoundLoc0, eBoundLoc0) = unitScale * unitScale;
    covariance(eBoundLoc0, eBoundDir0) = 0.9 * unitScale;
    covariance(eBoundDir0, eBoundLoc0) = 0.9 * unitScale;
    const BoundMatrix original = covariance;
    const auto result = detail::updateMaterialCovariance(
        covariance, noise, 0., NoiseUpdateMode::removeNoise);
    BOOST_REQUIRE(!result.ok());
    BOOST_CHECK(result.error() ==
                MaterialInteractionError::InvalidCovarianceAfterNoiseRemoval);
    BOOST_CHECK_EQUAL(covariance, original);
  }

  BoundMatrix covariance = BoundMatrix::Zero();
  noise(0, 0) = 1e-30;
  BOOST_CHECK(!detail::updateMaterialCovariance(covariance, noise, 0.,
                                                NoiseUpdateMode::removeNoise)
                   .ok());
  BOOST_CHECK(covariance.isZero(0.));
  BOOST_CHECK(!detail::updateMaterialCovariance(covariance,
                                                Matrix<2, 2>::Zero(), 1e-30,
                                                NoiseUpdateMode::removeNoise)
                   .ok());
  BOOST_CHECK(covariance.isZero(0.));
}

BOOST_AUTO_TEST_CASE(NoiseRemovalWithWidelySeparatedScales) {
  for (double unitScale : {1e-150, 1e-12, 1., 1e12, 1e150}) {
    BoundMatrix expected = BoundMatrix::Identity();
    expected(eBoundLoc0, eBoundLoc0) = unitScale * unitScale;
    expected(eBoundLoc0, eBoundDir0) = 0.5 * unitScale;
    expected(eBoundDir0, eBoundLoc0) = 0.5 * unitScale;
    BoundMatrix covariance = expected;
    Matrix<2, 2> noise;
    noise << 0.25, -0.125, -0.125, 0.25;
    covariance.block<2, 2>(eBoundDir0, eBoundDir0) += noise;
    covariance(eBoundQOverP, eBoundQOverP) += 0.125;
    const auto result = detail::updateMaterialCovariance(
        covariance, noise, 0.125, NoiseUpdateMode::removeNoise);
    BOOST_REQUIRE(result.ok());
    // Entrywise comparison avoids a large coordinate scale hiding errors.
    BOOST_CHECK_EQUAL(covariance, expected);
  }
}

}  // namespace ActsTests
