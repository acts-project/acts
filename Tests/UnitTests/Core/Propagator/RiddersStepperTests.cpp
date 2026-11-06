// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/RiddersStepper.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <numbers>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using Stepper = Experimental::RiddersStepper<EigenStepper<>>;

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

// A transport anchors the state, so a second transport at the same point does
// not change the covariance again
BOOST_AUTO_TEST_CASE(ridders_stepper_repeated_transport) {
  const GeometryContext geoCtx = GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext magCtx;

  Stepper stepper(std::make_shared<ConstantBField>(Vector3(0, 0, 2_T)));

  BoundMatrix cov = BoundMatrix::Identity();
  cov(eBoundQOverP, eBoundQOverP) = 1e-4;
  const auto start = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), 0.3, 1.2, 1. / 1_GeV, cov, ParticleHypothesis::pion());

  Stepper::Options options(geoCtx, magCtx);
  options.maxStepSize = 10_mm;
  auto state = stepper.makeState(options);
  stepper.initialize(state, start);

  for (int i = 0; i < 20; ++i) {
    BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
  }

  const auto surface = CurvilinearSurface(stepper.position(state) +
                                              5_mm * stepper.direction(state),
                                          stepper.direction(state))
                           .planeSurface();

  // Step all states onto the surface, like the navigator does
  bool onSurface = false;
  for (int i = 0; i < 10 && !onSurface; ++i) {
    onSurface =
        stepper.updateSurfaceStatus(
            state, *surface, 0, Direction::Forward(),
            BoundaryTolerance::Infinite(), s_onSurfaceTolerance,
            ConstrainedStep::Type::Navigator) == IntersectionStatus::onSurface;
    if (!onSurface) {
      BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
    }
  }
  BOOST_REQUIRE(onSurface);

  const auto firstJacobian = stepper.transportToBound(state, *surface);
  BOOST_REQUIRE(firstJacobian.ok());
  const BoundMatrix firstCovariance = stepper.covariance(state).value();
  BOOST_CHECK(!firstJacobian->isIdentity(1e-3));

  const auto secondJacobian = stepper.transportToBound(state, *surface);
  BOOST_REQUIRE(secondJacobian.ok());
  CHECK_CLOSE_ABS(*secondJacobian, BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_COVARIANCE(stepper.covariance(state).value(), firstCovariance,
                         1e-6);
}

BOOST_AUTO_TEST_CASE(ridders_matches_analytic_across_phi_seam) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  const MagneticFieldContext mctx;
  const auto field = std::make_shared<ConstantBField>(Vector3::Zero());
  Stepper::Config config;
  config.stepperConfig.bField = field;
  config.parameterVariation =
      std::make_shared<Experimental::DeltaBoundParameterVariationGenerator>(
          std::vector<double>{-1e-3, 1e-3});
  Stepper numerical(config);
  EigenStepper<> analytic(field);
  Stepper::Options options(gctx, mctx);
  options.maxStepSize = 10_mm;
  const auto source = Surface::makeShared<PlaneSurface>(Transform3::Identity());

  for (double phi : {-std::numbers::pi + 5e-4, std::numbers::pi - 5e-4}) {
    constexpr double theta = 1.;
    BOOST_TEST_CONTEXT("phi=" << phi << ", theta=" << theta) {
      BoundVector parameters = BoundVector::Zero();
      parameters[eBoundDir0] = phi;
      parameters[eBoundDir1] = theta;
      parameters[eBoundQOverP] = 1. / 1_GeV;
      BoundMatrix covariance = BoundMatrix::Identity();
      covariance(eBoundQOverP, eBoundQOverP) = 1e-4;
      covariance(eBoundLoc0, eBoundDir1) = .2;
      covariance(eBoundDir1, eBoundLoc0) = .2;
      const BoundTrackParameters start(source, parameters, covariance,
                                       ParticleHypothesis::pion());
      auto numericalState = numerical.makeState(options);
      auto analyticState = analytic.makeState(options);
      numerical.initialize(numericalState, start);
      analytic.initialize(analyticState, start);

      // Propagate to a parallel plane; secondary rays have different path
      // lengths. This exercises the surface-intersection derivatives too.
      Transform3 placement = Transform3::Identity();
      placement.translation() = 10_mm * start.direction(gctx);
      const auto target = Surface::makeShared<PlaneSurface>(placement);
      const auto reachSurface = [&](auto& stepper, auto& state) {
        for (int i = 0; i < 10; ++i) {
          const auto status = stepper.updateSurfaceStatus(
              state, *target, 0, Direction::Forward(),
              BoundaryTolerance::Infinite(), 1e-12,
              ConstrainedStep::Type::Navigator);
          if (status == IntersectionStatus::onSurface) {
            return true;
          }
          if (!stepper.step(state, Direction::Forward(), nullptr).ok()) {
            return false;
          }
        }
        return false;
      };
      BOOST_REQUIRE(reachSurface(numerical, numericalState));
      BOOST_REQUIRE(reachSurface(analytic, analyticState));
      const auto numericJac =
          numerical.transportToBound(numericalState, *target);
      const auto analyticJac =
          analytic.transportToBound(analyticState, *target);
      BOOST_REQUIRE(numericJac.ok());
      BOOST_REQUIRE(analyticJac.ok());
      CHECK_CLOSE_ABS(*numericJac, *analyticJac, 1e-4);
      CHECK_CLOSE_ABS(*numerical.covariance(numericalState),
                      *analytic.covariance(analyticState), 5e-3);
      const auto repeated = numerical.transportToBound(numericalState, *target);
      BOOST_REQUIRE(repeated.ok());
      CHECK_CLOSE_ABS(*repeated, BoundMatrix::Identity(), 1e-6);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
