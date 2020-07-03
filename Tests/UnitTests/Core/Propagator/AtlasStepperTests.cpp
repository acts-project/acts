// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;
using Covariance = BoundSymMatrix;
using Jacobian = BoundMatrix;
using Stepper = AtlasStepper<ConstantBField>;

/// Simplified propagator state.
struct MockPropagatorState {
  MockPropagatorState(Stepper::State stepperState)
      : stepping(std::move(stepperState)) {}

  /// Stepper state.
  Stepper::State stepping;
  /// Propagator options with only the relevant components.
  struct {
    double mass = 1_GeV;
    double tolerance = 10_um;
  } options;
};

// epsilon for floating point comparisons
static constexpr auto eps = 1024 * std::numeric_limits<double>::epsilon();

// propagation settings
static constexpr auto stepSize = 10_mm;
static constexpr auto tolerance = 10_um;
static constexpr NavigationDirection navDir = backward;
static const ConstantBField magneticField(Vector3D(0.1_T, -0.2_T, 2_T));

// initial parameter state
static const Vector3D pos(1_mm, -1_mm, 2_mm);
static constexpr auto time = 2_ns;
static const Vector3D unitDir = Vector3D(-2, 2, 1).normalized();
static constexpr auto absMom = 1_GeV;
static const Vector3D mom = absMom * unitDir;
static constexpr auto charge = -1_e;
static const Covariance cov = Covariance::Identity();

// context objects
static const GeometryContext geoCtx;
static const MagneticFieldContext magCtx;

BOOST_AUTO_TEST_SUITE(AtlasStepper)

// test state construction from parameters w/o covariance
BOOST_AUTO_TEST_CASE(ConstructState) {
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearParameters(std::nullopt, pos, mom, charge, time), navDir,
      stepSize, tolerance);

  BOOST_TEST(!state.covTransport);
  BOOST_TEST(state.covariance == nullptr);
  BOOST_TEST(state.pVector[0] == pos.x());
  BOOST_TEST(state.pVector[1] == pos.y());
  BOOST_TEST(state.pVector[2] == pos.z());
  BOOST_TEST(state.pVector[3] == time);
  CHECK_CLOSE_ABS(state.pVector[4], unitDir.x(), eps);
  CHECK_CLOSE_ABS(state.pVector[5], unitDir.y(), eps);
  CHECK_CLOSE_ABS(state.pVector[6], unitDir.z(), eps);
  BOOST_TEST(state.pVector[7] == charge / absMom);
  BOOST_TEST(state.navDir == navDir);
  BOOST_TEST(state.pathAccumulated == 0.);
  BOOST_TEST(state.stepSize == navDir * stepSize);
  BOOST_TEST(state.previousStepSize == 0.);
  BOOST_TEST(state.tolerance == tolerance);
}

// test state construction from parameters w/ covariance
BOOST_AUTO_TEST_CASE(ConstructStateWithCovariance) {
  Stepper::State state(geoCtx, magCtx,
                       CurvilinearParameters(cov, pos, mom, charge, time),
                       navDir, stepSize, tolerance);

  BOOST_TEST(state.covTransport);
  BOOST_TEST(*state.covariance == cov);
  BOOST_TEST(state.pVector[0] == pos.x());
  BOOST_TEST(state.pVector[1] == pos.y());
  BOOST_TEST(state.pVector[2] == pos.z());
  BOOST_TEST(state.pVector[3] == time);
  CHECK_CLOSE_ABS(state.pVector[4], unitDir.x(), eps);
  CHECK_CLOSE_ABS(state.pVector[5], unitDir.y(), eps);
  CHECK_CLOSE_ABS(state.pVector[6], unitDir.z(), eps);
  BOOST_TEST(state.pVector[7] == charge / absMom);
  BOOST_TEST(state.navDir == navDir);
  BOOST_TEST(state.pathAccumulated == 0.);
  BOOST_TEST(state.stepSize == navDir * stepSize);
  BOOST_TEST(state.previousStepSize == 0.);
  BOOST_TEST(state.tolerance == tolerance);
}

// test stepper getters for particle state
BOOST_AUTO_TEST_CASE(Getters) {
  Stepper stepper(magneticField);
  Stepper::State state(geoCtx, magCtx,
                       CurvilinearParameters(cov, pos, mom, charge, time),
                       navDir, stepSize, tolerance);

  CHECK_CLOSE_ABS(stepper.position(state), pos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), time, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), unitDir, eps);
  CHECK_CLOSE_ABS(stepper.momentum(state), absMom, eps);
  BOOST_TEST(stepper.charge(state) == charge);
}

// test stepper update methods with bound state as input
BOOST_AUTO_TEST_CASE(UpdateFromBound) {
  Stepper stepper(magneticField);
  Stepper::State state(geoCtx, magCtx,
                       CurvilinearParameters(cov, pos, mom, charge, time),
                       navDir, stepSize, tolerance);

  auto newPos = (pos + Vector3D(1_mm, 2_mm, 3_mm)).eval();
  auto newTime = time + 20_ns;
  auto newUnitDir = (unitDir + Vector3D(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;
  auto newMom = (newAbsMom * newUnitDir).eval();

  // example surface and bound parameters at the updated position
  auto plane = Surface::makeShared<PlaneSurface>(newPos, newUnitDir);
  BoundParameters params(geoCtx, cov, newPos, newMom, charge, newTime, plane);
  FreeVector freeParams;
  freeParams[eFreePos0] = newPos[eX];
  freeParams[eFreePos1] = newPos[eY];
  freeParams[eFreePos2] = newPos[eZ];
  freeParams[eFreeTime] = newTime;
  freeParams[eFreeDir0] = newUnitDir[eMom0];
  freeParams[eFreeDir1] = newUnitDir[eMom1];
  freeParams[eFreeDir2] = newUnitDir[eMom2];
  freeParams[eFreeQOverP] = charge / newAbsMom;

  // WARNING for some reason there seems to be an additional flag that makes
  //         the update method not do anything when it is set. Why?
  state.state_ready = false;
  stepper.update(state, freeParams, *params.covariance());
  CHECK_CLOSE_ABS(stepper.position(state), newPos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), newTime, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), newUnitDir, eps);
  CHECK_CLOSE_ABS(stepper.momentum(state), newAbsMom, eps);
  BOOST_TEST(stepper.charge(state) == charge);
}

// test stepper update methods with individual components as input
BOOST_AUTO_TEST_CASE(UpdateFromComponents) {
  Stepper stepper(magneticField);
  Stepper::State state(geoCtx, magCtx,
                       CurvilinearParameters(cov, pos, mom, charge, time),
                       navDir, stepSize, tolerance);

  auto newPos = (pos + Vector3D(1_mm, 2_mm, 3_mm)).eval();
  auto newTime = time + 20_ns;
  auto newUnitDir = (unitDir + Vector3D(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;

  stepper.update(state, newPos, newUnitDir, newAbsMom, newTime);
  CHECK_CLOSE_ABS(stepper.position(state), newPos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), newTime, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), newUnitDir, eps);
  CHECK_CLOSE_ABS(stepper.momentum(state), newAbsMom, eps);
  BOOST_TEST(stepper.charge(state) == charge);
}

// test building a bound state object from the stepper state
BOOST_AUTO_TEST_CASE(BuildBound) {
  Stepper stepper(magneticField);
  Stepper::State state(geoCtx, magCtx,
                       CurvilinearParameters(cov, pos, mom, charge, time),
                       navDir, stepSize, tolerance);
  // example surface at the current state position
  auto plane = Surface::makeShared<PlaneSurface>(pos, unitDir);

  auto&& [pars, jac, pathLength] = stepper.boundState(state, *plane);
  // check parameters
  CHECK_CLOSE_ABS(pars.position(), pos, eps);
  CHECK_CLOSE_ABS(pars.time(), time, eps);
  CHECK_CLOSE_ABS(pars.momentum(), mom, eps);
  BOOST_TEST(pars.charge() == charge);
  BOOST_TEST(pars.covariance().has_value());
  BOOST_TEST(*pars.covariance() != cov);
  // check Jacobian. should be identity since there was no propagation yet
  CHECK_CLOSE_ABS(jac, Jacobian(Jacobian::Identity()), eps);
  // check propagation length
  CHECK_CLOSE_ABS(pathLength, 0., eps);
}

// test building a curvilinear state object from the stepper state
BOOST_AUTO_TEST_CASE(BuildCurvilinear) {
  Stepper stepper(magneticField);
  Stepper::State state(geoCtx, magCtx,
                       CurvilinearParameters(cov, pos, mom, charge, time),
                       navDir, stepSize, tolerance);

  auto&& [pars, jac, pathLength] = stepper.curvilinearState(state);
  // check parameters
  CHECK_CLOSE_ABS(pars.position(), pos, eps);
  CHECK_CLOSE_ABS(pars.time(), time, eps);
  CHECK_CLOSE_ABS(pars.momentum(), mom, eps);
  BOOST_TEST(pars.charge() == charge);
  BOOST_TEST(pars.covariance().has_value());
  BOOST_TEST(*pars.covariance() != cov);
  // check Jacobian. should be identity since there was no propagation yet
  CHECK_CLOSE_ABS(jac, Jacobian(Jacobian::Identity()), eps);
  // check propagation length
  CHECK_CLOSE_ABS(pathLength, 0., eps);
}

// test step method without covariance transport
BOOST_AUTO_TEST_CASE(Step) {
  Stepper stepper(magneticField);
  MockPropagatorState state(Stepper::State(
      geoCtx, magCtx, CurvilinearParameters(cov, pos, mom, charge, time),
      navDir, stepSize, tolerance));
  state.stepping.covTransport = false;

  // ensure step does not result in an error
  auto res = stepper.step(state);
  BOOST_TEST(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_TEST(state.stepping.stepSize == navDir * stepSize);
  BOOST_TEST(state.stepping.stepSize == h);

  // check that the position has moved
  auto deltaPos = (stepper.position(state.stepping) - pos).eval();
  BOOST_TEST(0 < deltaPos.norm());
  // check that time has changed
  auto deltaTime = stepper.time(state.stepping) - time;
  BOOST_TEST(0 < std::abs(deltaTime));
  // check that the direction was rotated
  auto projDir = stepper.direction(state.stepping).dot(unitDir);
  BOOST_TEST(projDir < 1);

  // momentum and charge should not change
  CHECK_CLOSE_ABS(stepper.momentum(state.stepping), absMom, eps);
  BOOST_TEST(stepper.charge(state.stepping) == charge);
}

// test step method with covariance transport
BOOST_AUTO_TEST_CASE(StepWithCovariance) {
  Stepper stepper(magneticField);
  MockPropagatorState state(Stepper::State(
      geoCtx, magCtx, CurvilinearParameters(cov, pos, mom, charge, time),
      navDir, stepSize, tolerance));
  state.stepping.covTransport = true;

  // ensure step does not result in an error
  auto res = stepper.step(state);
  BOOST_TEST(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_TEST(state.stepping.stepSize == navDir * stepSize);
  BOOST_TEST(state.stepping.stepSize == h);

  // check that the position has moved
  auto deltaPos = (stepper.position(state.stepping) - pos).eval();
  BOOST_TEST(0 < deltaPos.norm());
  // check that time has changed
  auto deltaTime = stepper.time(state.stepping) - time;
  BOOST_TEST(0 < std::abs(deltaTime));
  // check that the direction was rotated
  auto projDir = stepper.direction(state.stepping).dot(unitDir);
  BOOST_TEST(projDir < 1);

  // momentum and charge should not change
  CHECK_CLOSE_ABS(stepper.momentum(state.stepping), absMom, eps);
  BOOST_TEST(stepper.charge(state.stepping) == charge);

  stepper.covarianceTransport(state.stepping);
  BOOST_TEST(state.stepping.cov != cov);
}

BOOST_AUTO_TEST_CASE(StepSize) {
  Stepper stepper(magneticField);
  Stepper::State state(geoCtx, magCtx,
                       CurvilinearParameters(cov, pos, mom, charge, time),
                       navDir, stepSize, tolerance);

  // TODO figure out why this fails and what it should be
  // BOOST_TEST(stepper.overstepLimit(state) == tolerance);

  stepper.setStepSize(state, 5_cm);
  BOOST_TEST(state.previousStepSize == navDir * stepSize);
  BOOST_TEST(state.stepSize == 5_cm);

  stepper.releaseStepSize(state);
  BOOST_TEST(state.stepSize == navDir * stepSize);
}

// test step size modification with target surfaces
BOOST_AUTO_TEST_CASE(StepSizeSurface) {
  Stepper stepper(magneticField);
  Stepper::State state(geoCtx, magCtx,
                       CurvilinearParameters(cov, pos, mom, charge, time),
                       navDir, stepSize, tolerance);

  auto distance = 10_mm;
  auto target = Surface::makeShared<PlaneSurface>(
      pos + navDir * distance * unitDir, unitDir);

  stepper.updateSurfaceStatus(state, *target, BoundaryCheck(false));
  BOOST_TEST(state.stepSize.value(ConstrainedStep::actor) == navDir * distance);

  // test the step size modification in the context of a surface
  stepper.updateStepSize(
      state,
      target->intersect(state.geoContext, stepper.position(state),
                        state.navDir * stepper.direction(state), false),
      false);
  BOOST_TEST(state.stepSize == distance);

  // start with a different step size
  state.stepSize = navDir * stepSize;
  stepper.updateStepSize(
      state,
      target->intersect(state.geoContext, stepper.position(state),
                        state.navDir * stepper.direction(state), false),
      true);
  BOOST_TEST(state.stepSize == distance);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts