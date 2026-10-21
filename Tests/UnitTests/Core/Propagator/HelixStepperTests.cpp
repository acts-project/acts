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
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/HelixStepper.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

namespace ActsTests {

namespace {

const auto tgContext = GeometryContext::dangerouslyDefaultConstruct();
const MagneticFieldContext mfContext;

/// Field with a constant gradient of Bz along x, to force step rejections.
class GradientField final : public MagneticFieldProvider {
 public:
  struct Cache {
    explicit Cache(const MagneticFieldContext& /*mctx*/) {}
  };

  GradientField(double bz, double dBzdx) : m_bz(bz), m_dBzdx(dBzdx) {}

  MagneticFieldProvider::Cache makeCache(
      const MagneticFieldContext& mctx) const override {
    return MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Result<Vector3> getField(
      const Vector3& position,
      MagneticFieldProvider::Cache& /*cache*/) const override {
    return Result<Vector3>::success(
        Vector3(0., 0., m_bz + m_dBzdx * position.x()));
  }

 private:
  double m_bz;
  double m_dBzdx;
};

BoundTrackParameters makeParameters(const Vector3& pos, const Vector3& dir,
                                    double qop,
                                    std::optional<BoundMatrix> cov) {
  return BoundTrackParameters::createCurvilinear(
      makeVector4(pos, 0.), dir, qop, cov, ParticleHypothesis::pion());
}

/// Run one step of fixed size from the free parameters and return them.
FreeVector stepFrom(const HelixStepper& stepper, const FreeVector& start,
                    double h, FreeMatrix* jacTransport = nullptr) {
  HelixStepper::Options options(tgContext, mfContext);
  options.initialStepSize = std::abs(h);
  options.maxStepSize = std::abs(h);
  HelixStepper::State state = stepper.makeState(options);
  stepper.initialize(state, makeParameters(Vector3::Zero(), Vector3::UnitX(),
                                           1., BoundMatrix::Identity()));
  state.pars = start;
  state.field.reset();
  auto res = stepper.step(state, Direction::fromScalar(h), nullptr);
  BOOST_REQUIRE(res.ok());
  CHECK_CLOSE_REL(*res, h, 1e-15);
  if (jacTransport != nullptr) {
    *jacTransport = state.jacTransport;
  }
  return state.pars;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

BOOST_AUTO_TEST_CASE(helix_stepper_state_test) {
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 2_T));
  HelixStepper stepper(bField);

  HelixStepper::Options options(tgContext, mfContext);
  options.maxStepSize = 123.;
  HelixStepper::State state = stepper.makeState(options);

  Vector3 pos(1., 2., 3.);
  Vector3 dir = Vector3(4., 5., 6.).normalized();
  stepper.initialize(state,
                     makeParameters(pos, dir, -1. / 8_GeV, std::nullopt));

  BOOST_CHECK(!stepper.hasCovariance(state));
  BOOST_CHECK(!state.field.has_value());
  CHECK_CLOSE_ABS(stepper.position(state), pos, 1e-12);
  CHECK_CLOSE_ABS(stepper.direction(state), dir, 1e-12);
  CHECK_CLOSE_REL(stepper.absoluteMomentum(state), 8_GeV, 1e-12);
  BOOST_CHECK_EQUAL(stepper.charge(state), -1.);
  BOOST_CHECK_EQUAL(state.stepSize.value(), 123.);

  stepper.initialize(
      state, makeParameters(pos, dir, -1. / 8_GeV, BoundMatrix::Identity()));
  BOOST_CHECK(stepper.hasCovariance(state));
  BOOST_CHECK_NE(state.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(state.jacTransport, FreeMatrix::Identity());
}

/// A full turn of the helix returns to the start, shifted along the field.
BOOST_AUTO_TEST_CASE(helix_stepper_full_turn) {
  const double bz = 2_T;
  const double p = 1_GeV;
  const double qop = 1. / p;
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., bz));
  HelixStepper stepper(bField);

  const Vector3 pos(10., -20., 30.);
  const Vector3 dir = Vector3(1., 1., 0.5).normalized();
  // The direction turns about the field once per this path length
  const double turn = 2. * std::numbers::pi / (std::abs(qop) * bz);

  for (int nSteps : {1, 7}) {
    HelixStepper::Options options(tgContext, mfContext);
    options.initialStepSize = turn / nSteps;
    options.maxStepSize = turn / nSteps;
    HelixStepper::State state = stepper.makeState(options);
    stepper.initialize(state, makeParameters(pos, dir, qop, std::nullopt));
    for (int i = 0; i < nSteps; ++i) {
      BOOST_REQUIRE(stepper.step(state, Direction::Forward(), nullptr).ok());
    }
    CHECK_CLOSE_ABS(stepper.position(state),
                    Vector3(pos + Vector3(0., 0., turn * dir.z())), 1e-9);
    CHECK_CLOSE_ABS(stepper.direction(state), dir, 1e-12);
    BOOST_CHECK_EQUAL(state.statistics.nRejectedSteps, 0u);
  }
}

/// The helix with a vanishing field is a straight line.
BOOST_AUTO_TEST_CASE(helix_stepper_straight_limit) {
  auto bField = std::make_shared<ConstantBField>(Vector3::Zero());
  HelixStepper helix(bField);
  StraightLineStepper line;

  const double h = 123.;
  HelixStepper::Options helixOptions(tgContext, mfContext);
  helixOptions.maxStepSize = h;
  StraightLineStepper::Options lineOptions(tgContext, mfContext);
  lineOptions.maxStepSize = h;

  auto params =
      makeParameters(Vector3(1., 2., 3.), Vector3(4., 5., 6.).normalized(),
                     -1. / 2_GeV, BoundMatrix::Identity());
  HelixStepper::State helixState = helix.makeState(helixOptions);
  helix.initialize(helixState, params);
  StraightLineStepper::State lineState = line.makeState(lineOptions);
  line.initialize(lineState, params);

  BOOST_REQUIRE(helix.step(helixState, Direction::Forward(), nullptr).ok());
  BOOST_REQUIRE(line.step(lineState, Direction::Forward(), nullptr).ok());

  CHECK_CLOSE_ABS(helixState.pars, lineState.pars, 1e-12);
  CHECK_CLOSE_ABS(helixState.jacTransport, lineState.jacTransport, 1e-12);
  CHECK_CLOSE_ABS(helixState.derivative, lineState.derivative, 1e-12);
}

/// The analytic step jacobian matches central finite differences.
BOOST_AUTO_TEST_CASE(helix_stepper_jacobian) {
  auto bField = std::make_shared<ConstantBField>(Vector3(0.1_T, -0.2_T, 2_T));
  HelixStepper stepper(bField);

  FreeVector start = FreeVector::Zero();
  start.segment<3>(eFreePos0) = Vector3(1., 2., 3.);
  start[eFreeTime] = 4.;
  start.segment<3>(eFreeDir0) = Vector3(4., -5., 6.).normalized();
  start[eFreeQOverP] = -1. / 0.7_GeV;

  // Small and large turning angles, forward and backward
  for (double h : {1e-3, 50., 1500., -400.}) {
    FreeMatrix jac;
    stepFrom(stepper, start, h, &jac);

    std::array<double, eFreeSize> deltas{};
    deltas.fill(1e-6);
    deltas[eFreeQOverP] = 1e-6 * std::abs(start[eFreeQOverP]);

    FreeMatrix numeric;
    for (std::size_t j = 0; j < eFreeSize; ++j) {
      FreeVector plus = start;
      FreeVector minus = start;
      plus[j] += deltas[j];
      minus[j] -= deltas[j];
      numeric.col(j) =
          (stepFrom(stepper, plus, h) - stepFrom(stepper, minus, h)) /
          (2. * deltas[j]);
    }

    BOOST_TEST_CONTEXT("h = " << h) {
      CHECK_CLOSE_OR_SMALL(jac, numeric, 1e-5, 1e-6);
    }
  }
}

/// A field gradient shortens the step.
BOOST_AUTO_TEST_CASE(helix_stepper_step_size_control) {
  const auto params = makeParameters(Vector3::Zero(), Vector3::UnitX(),
                                     1. / 1_GeV, std::nullopt);

  HelixStepper::Options options(tgContext, mfContext);
  options.stepTolerance = 1e-4;
  options.initialStepSize = 1_m;
  options.maxStepSize = 1_m;

  HelixStepper gradient(std::make_shared<GradientField>(2_T, 1_T / 1_m));
  HelixStepper::State gradientState = gradient.makeState(options);
  gradient.initialize(gradientState, params);
  auto gradientStep =
      gradient.step(gradientState, Direction::Forward(), nullptr);
  BOOST_REQUIRE(gradientStep.ok());
  BOOST_CHECK_LT(*gradientStep, 1_m);
  BOOST_CHECK_GT(gradientState.statistics.nRejectedSteps, 0u);
  // The field at the end of the step is cached for the next step
  BOOST_REQUIRE(gradientState.field.has_value());
  CHECK_CLOSE_ABS(
      *gradientState.field,
      Vector3(0., 0., 2_T + 1_T / 1_m * gradient.position(gradientState).x()),
      1e-15);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
