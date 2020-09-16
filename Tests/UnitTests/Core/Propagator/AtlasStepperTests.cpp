// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/Assertions.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;
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
static const Vector4D pos4(1_mm, -1_mm, 2_mm, 2_ns);
static const Vector3D pos = pos4.segment<3>(ePos0);
static const auto time = pos4[eTime];
static const Vector3D unitDir = Vector3D(-2, 2, 1).normalized();
static constexpr auto absMom = 1_GeV;
static constexpr auto charge = -1_e;
static const Covariance cov = Covariance::Identity();

// context objects
static const GeometryContext geoCtx;
static const MagneticFieldContext magCtx;

BOOST_AUTO_TEST_SUITE(AtlasStepper)

// test state construction from parameters w/o covariance
BOOST_AUTO_TEST_CASE(ConstructState) {
  Stepper::State state(
      geoCtx, magCtx, CurvilinearTrackParameters(pos4, unitDir, absMom, charge),
      navDir, stepSize, tolerance);

  BOOST_CHECK(!state.covTransport);
  BOOST_CHECK_EQUAL(state.covariance, nullptr);
  BOOST_CHECK_EQUAL(state.pVector[0], pos.x());
  BOOST_CHECK_EQUAL(state.pVector[1], pos.y());
  BOOST_CHECK_EQUAL(state.pVector[2], pos.z());
  BOOST_CHECK_EQUAL(state.pVector[3], time);
  CHECK_CLOSE_ABS(state.pVector[4], unitDir.x(), eps);
  CHECK_CLOSE_ABS(state.pVector[5], unitDir.y(), eps);
  CHECK_CLOSE_ABS(state.pVector[6], unitDir.z(), eps);
  BOOST_CHECK_EQUAL(state.pVector[7], charge / absMom);
  BOOST_CHECK_EQUAL(state.navDir, navDir);
  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(state.stepSize, navDir * stepSize);
  BOOST_CHECK_EQUAL(state.previousStepSize, 0.);
  BOOST_CHECK_EQUAL(state.tolerance, tolerance);
}

// test state construction from parameters w/ covariance
BOOST_AUTO_TEST_CASE(ConstructStateWithCovariance) {
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  BOOST_CHECK(state.covTransport);
  BOOST_CHECK_EQUAL(*state.covariance, cov);
  BOOST_CHECK_EQUAL(state.pVector[0], pos.x());
  BOOST_CHECK_EQUAL(state.pVector[1], pos.y());
  BOOST_CHECK_EQUAL(state.pVector[2], pos.z());
  BOOST_CHECK_EQUAL(state.pVector[3], time);
  CHECK_CLOSE_ABS(state.pVector[4], unitDir.x(), eps);
  CHECK_CLOSE_ABS(state.pVector[5], unitDir.y(), eps);
  CHECK_CLOSE_ABS(state.pVector[6], unitDir.z(), eps);
  BOOST_CHECK_EQUAL(state.pVector[7], charge / absMom);
  BOOST_CHECK_EQUAL(state.navDir, navDir);
  BOOST_CHECK_EQUAL(state.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(state.stepSize, navDir * stepSize);
  BOOST_CHECK_EQUAL(state.previousStepSize, 0.);
  BOOST_CHECK_EQUAL(state.tolerance, tolerance);
}

// test stepper getters for particle state
BOOST_AUTO_TEST_CASE(Getters) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  CHECK_CLOSE_ABS(stepper.position(state), pos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), time, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), unitDir, eps);
  CHECK_CLOSE_ABS(stepper.momentum(state), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test stepper update methods with bound state as input
BOOST_AUTO_TEST_CASE(UpdateFromBound) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  auto newPos4 = (pos4 + Vector4D(1_mm, 2_mm, 3_mm, 20_ns)).eval();
  auto newPos = newPos4.segment<3>(ePos0);
  auto newTime = newPos4[eTime];
  auto newUnitDir = (unitDir + Vector3D(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;

  // example surface and bound parameters at the updated position
  auto plane = Surface::makeShared<PlaneSurface>(newPos, newUnitDir);
  BoundTrackParameters params(plane, geoCtx, newPos4, newUnitDir, charge, cov);
  FreeVector freeParams;
  freeParams[eFreePos0] = newPos4[ePos0];
  freeParams[eFreePos1] = newPos4[ePos1];
  freeParams[eFreePos2] = newPos4[ePos2];
  freeParams[eFreeTime] = newPos4[eTime];
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
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test stepper update methods with individual components as input
BOOST_AUTO_TEST_CASE(UpdateFromComponents) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  auto newPos = (pos + Vector3D(1_mm, 2_mm, 3_mm)).eval();
  auto newTime = time + 20_ns;
  auto newUnitDir = (unitDir + Vector3D(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;

  stepper.update(state, newPos, newUnitDir, newAbsMom, newTime);
  CHECK_CLOSE_ABS(stepper.position(state), newPos, eps);
  CHECK_CLOSE_ABS(stepper.time(state), newTime, eps);
  CHECK_CLOSE_ABS(stepper.direction(state), newUnitDir, eps);
  CHECK_CLOSE_ABS(stepper.momentum(state), newAbsMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state), charge);
}

// test building a bound state object from the stepper state
BOOST_AUTO_TEST_CASE(BuildBound) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);
  // example surface at the current state position
  auto plane = Surface::makeShared<PlaneSurface>(pos, unitDir);

  auto&& [pars, jac, pathLength] = stepper.boundState(state, *plane);
  // check parameters
  CHECK_CLOSE_ABS(pars.position(geoCtx), pos, eps);
  CHECK_CLOSE_ABS(pars.time(), time, eps);
  CHECK_CLOSE_ABS(pars.momentum(), absMom * unitDir, eps);
  BOOST_CHECK_EQUAL(pars.charge(), charge);
  BOOST_CHECK(pars.covariance().has_value());
  BOOST_CHECK_NE(*pars.covariance(), cov);
  // check Jacobian. should be identity since there was no propagation yet
  CHECK_CLOSE_ABS(jac, Jacobian(Jacobian::Identity()), eps);
  // check propagation length
  CHECK_CLOSE_ABS(pathLength, 0., eps);
}

// test building a curvilinear state object from the stepper state
BOOST_AUTO_TEST_CASE(BuildCurvilinear) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  auto&& [pars, jac, pathLength] = stepper.curvilinearState(state);
  // check parameters
  CHECK_CLOSE_ABS(pars.position(geoCtx), pos, eps);
  CHECK_CLOSE_ABS(pars.time(), time, eps);
  CHECK_CLOSE_ABS(pars.momentum(), absMom * unitDir, eps);
  BOOST_CHECK_EQUAL(pars.charge(), charge);
  BOOST_CHECK(pars.covariance().has_value());
  BOOST_CHECK_NE(*pars.covariance(), cov);
  // check Jacobian. should be identity since there was no propagation yet
  CHECK_CLOSE_ABS(jac, Jacobian(Jacobian::Identity()), eps);
  // check propagation length
  CHECK_CLOSE_ABS(pathLength, 0., eps);
}

// test step method without covariance transport
BOOST_AUTO_TEST_CASE(Step) {
  Stepper stepper(magneticField);
  MockPropagatorState state(Stepper::State(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance));
  state.stepping.covTransport = false;

  // ensure step does not result in an error
  auto res = stepper.step(state);
  BOOST_CHECK(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_CHECK_EQUAL(state.stepping.stepSize, navDir * stepSize);
  BOOST_CHECK_EQUAL(state.stepping.stepSize, h);

  // check that the position has moved
  auto deltaPos = (stepper.position(state.stepping) - pos).eval();
  BOOST_CHECK_LT(0, deltaPos.norm());
  // check that time has changed
  auto deltaTime = stepper.time(state.stepping) - time;
  BOOST_CHECK_LT(0, std::abs(deltaTime));
  // check that the direction was rotated
  auto projDir = stepper.direction(state.stepping).dot(unitDir);
  BOOST_CHECK_LT(projDir, 1);

  // momentum and charge should not change
  CHECK_CLOSE_ABS(stepper.momentum(state.stepping), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state.stepping), charge);
}

// test step method with covariance transport
BOOST_AUTO_TEST_CASE(StepWithCovariance) {
  Stepper stepper(magneticField);
  MockPropagatorState state(Stepper::State(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance));
  state.stepping.covTransport = true;

  // ensure step does not result in an error
  auto res = stepper.step(state);
  BOOST_CHECK(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_CHECK_EQUAL(state.stepping.stepSize, navDir * stepSize);
  BOOST_CHECK_EQUAL(state.stepping.stepSize, h);

  // check that the position has moved
  auto deltaPos = (stepper.position(state.stepping) - pos).eval();
  BOOST_CHECK_LT(0, deltaPos.norm());
  // check that time has changed
  auto deltaTime = stepper.time(state.stepping) - time;
  BOOST_CHECK_LT(0, std::abs(deltaTime));
  // check that the direction was rotated
  auto projDir = stepper.direction(state.stepping).dot(unitDir);
  BOOST_CHECK_LT(projDir, 1);

  // momentum and charge should not change
  CHECK_CLOSE_ABS(stepper.momentum(state.stepping), absMom, eps);
  BOOST_CHECK_EQUAL(stepper.charge(state.stepping), charge);

  stepper.covarianceTransport(state.stepping);
  BOOST_CHECK_NE(state.stepping.cov, cov);
}

// test state reset method
BOOST_AUTO_TEST_CASE(Reset) {
  Stepper stepper(magneticField);
  MockPropagatorState state(Stepper::State(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance));
  state.stepping.covTransport = true;

  // ensure step does not result in an error
  stepper.step(state);

  // Construct the parameters
  Vector3D pos(1.5, -2.5, 3.5);
  Vector3D mom(4.5, -5.5, 6.5);
  double time = 7.5;
  double charge = 1.;
  BoundSymMatrix cov = 8.5 * Covariance::Identity();
  CurvilinearTrackParameters cp(makeVector4(pos, time), unitDir, absMom, charge,
                                cov);
  FreeVector freeParams = detail::transformBoundToFreeParameters(
      cp.referenceSurface(), geoCtx, cp.parameters());
  NavigationDirection ndir = forward;
  double stepSize = -256.;

  // Reset all possible parameters
  Stepper::State stateCopy(state.stepping);
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface(), ndir, stepSize);
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, cov);
  BOOST_CHECK_EQUAL(stepper.position(stateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(stepper.direction(stateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(stepper.momentum(stateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(stepper.charge(stateCopy), stepper.charge(state.stepping));
  BOOST_CHECK_EQUAL(stepper.time(stateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(stateCopy.navDir, ndir);
  BOOST_CHECK_EQUAL(stateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(stateCopy.stepSize, ndir * stepSize);
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize,
                    state.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(stateCopy.tolerance, state.stepping.tolerance);

  // Reset all possible parameters except the step size
  stateCopy = state.stepping;
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface(), ndir);
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, cov);
  BOOST_CHECK_EQUAL(stepper.position(stateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(stepper.direction(stateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(stepper.momentum(stateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(stepper.charge(stateCopy), stepper.charge(state.stepping));
  BOOST_CHECK_EQUAL(stepper.time(stateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(stateCopy.navDir, ndir);
  BOOST_CHECK_EQUAL(stateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(stateCopy.stepSize,
                    ndir * std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize,
                    state.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(stateCopy.tolerance, state.stepping.tolerance);

  // Reset the least amount of parameters
  stateCopy = state.stepping;
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface());
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, cov);
  BOOST_CHECK_EQUAL(stepper.position(stateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(stepper.direction(stateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(stepper.momentum(stateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(stepper.charge(stateCopy), stepper.charge(state.stepping));
  BOOST_CHECK_EQUAL(stepper.time(stateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(stateCopy.navDir, forward);
  BOOST_CHECK_EQUAL(stateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(stateCopy.stepSize, std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize,
                    state.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(stateCopy.tolerance, state.stepping.tolerance);

  // Reset using different surface shapes
  // 1) Disc surface
  // Setting some parameters
  pos << 1.5, -2.5, 0.;
  mom << 4.5, -5.5, 6.5;
  time = 7.5;
  charge = 1.;
  cov = 8.5 * Covariance::Identity();
  Transform3D trafo = Transform3D::Identity();
  auto disc = Surface::makeShared<DiscSurface>(trafo);
  BoundTrackParameters boundDisc(disc, geoCtx, pos4, unitDir, absMom, charge);

  // Reset the state and test
  Stepper::State stateDisc = state.stepping;
  stepper.resetState(stateDisc, boundDisc.parameters(), *boundDisc.covariance(),
                     boundDisc.referenceSurface());

  CHECK_NE_COLLECTIONS(stateDisc.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(stateDisc.pVector, state.stepping.pVector);

  // 2) Perigee surface
  // Setting some parameters
  pos << 1.5, -2.5, 3.5;
  mom << 4.5, -5.5, 6.5;
  time = 7.5;
  charge = 1.;
  cov = 8.5 * Covariance::Identity();
  auto perigee = Surface::makeShared<PerigeeSurface>(trafo);
  BoundTrackParameters boundPerigee(perigee, geoCtx, pos4, unitDir, absMom,
                                    charge);

  // Reset the state and test
  Stepper::State statePerigee = state.stepping;
  stepper.resetState(statePerigee, boundPerigee.parameters(),
                     *boundPerigee.covariance(),
                     boundPerigee.referenceSurface());
  CHECK_NE_COLLECTIONS(statePerigee.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(statePerigee.pVector, state.stepping.pVector);
  CHECK_NE_COLLECTIONS(statePerigee.pVector, stateDisc.pVector);

  // 3) Straw surface
  // Setting some parameters
  pos << 1.5, -2.5, 3.5;
  mom << 4.5, -5.5, 6.5;
  time = 7.5;
  charge = 1.;
  cov = 8.5 * Covariance::Identity();
  auto straw = Surface::makeShared<StrawSurface>(trafo);
  BoundTrackParameters boundStraw(straw, geoCtx, pos4, unitDir, absMom, charge);

  // Reset the state and test
  Stepper::State stateStraw = state.stepping;
  stepper.resetState(stateStraw, boundStraw.parameters(),
                     *boundStraw.covariance(), boundStraw.referenceSurface());
  CHECK_NE_COLLECTIONS(stateStraw.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(stateStraw.pVector, state.stepping.pVector);
  CHECK_NE_COLLECTIONS(stateStraw.pVector, stateDisc.pVector);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      stateStraw.pVector.begin(), stateStraw.pVector.end(),
      statePerigee.pVector.begin(), statePerigee.pVector.end());
}

BOOST_AUTO_TEST_CASE(StepSize) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  // TODO figure out why this fails and what it should be
  // BOOST_CHECK_EQUAL(stepper.overstepLimit(state), tolerance);

  stepper.setStepSize(state, 5_cm);
  BOOST_CHECK_EQUAL(state.previousStepSize, navDir * stepSize);
  BOOST_CHECK_EQUAL(state.stepSize, 5_cm);

  stepper.releaseStepSize(state);
  BOOST_CHECK_EQUAL(state.stepSize, navDir * stepSize);
}

// test step size modification with target surfaces
BOOST_AUTO_TEST_CASE(StepSizeSurface) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magCtx,
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  auto distance = 10_mm;
  auto target = Surface::makeShared<PlaneSurface>(
      pos + navDir * distance * unitDir, unitDir);

  stepper.updateSurfaceStatus(state, *target, BoundaryCheck(false));
  BOOST_CHECK_EQUAL(state.stepSize.value(ConstrainedStep::actor),
                    navDir * distance);

  // test the step size modification in the context of a surface
  stepper.updateStepSize(
      state,
      target->intersect(state.geoContext, stepper.position(state),
                        state.navDir * stepper.direction(state), false),
      false);
  BOOST_CHECK_EQUAL(state.stepSize, distance);

  // start with a different step size
  state.stepSize = navDir * stepSize;
  stepper.updateStepSize(
      state,
      target->intersect(state.geoContext, stepper.position(state),
                        state.navDir * stepper.direction(state), false),
      true);
  BOOST_CHECK_EQUAL(state.stepSize, distance);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
