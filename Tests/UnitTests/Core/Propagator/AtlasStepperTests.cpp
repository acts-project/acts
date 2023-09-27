// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/Assertions.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;
using Covariance = BoundSymMatrix;
using Jacobian = BoundMatrix;
using Stepper = Acts::AtlasStepper;

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

struct MockNavigator {};

static constexpr MockNavigator navigator;

// epsilon for floating point comparisons
static constexpr auto eps = 1024 * std::numeric_limits<double>::epsilon();

// propagation settings
static constexpr auto stepSize = 10_mm;
static constexpr auto tolerance = 10_um;
static constexpr NavigationDirection navDir = NavigationDirection::Backward;
static auto magneticField =
    std::make_shared<ConstantBField>(Vector3(0.1_T, -0.2_T, 2_T));

// initial parameter state
static const Vector4 pos4(1_mm, -1_mm, 2_mm, 2_ns);
static const Vector3 pos = pos4.segment<3>(ePos0);
static const auto time = pos4[eTime];
static const Vector3 unitDir = Vector3(-2, 2, 1).normalized();
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
      geoCtx, magneticField->makeCache(magCtx),
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge), navDir,
      stepSize, tolerance);

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
  BOOST_CHECK_EQUAL(state.stepSize.value(), navDir * stepSize);
  BOOST_CHECK_EQUAL(state.previousStepSize, 0.);
  BOOST_CHECK_EQUAL(state.tolerance, tolerance);
}

// test state construction from parameters w/ covariance
BOOST_AUTO_TEST_CASE(ConstructStateWithCovariance) {
  Stepper::State state(
      geoCtx, magneticField->makeCache(magCtx),
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
  BOOST_CHECK_EQUAL(state.stepSize.value(), navDir * stepSize);
  BOOST_CHECK_EQUAL(state.previousStepSize, 0.);
  BOOST_CHECK_EQUAL(state.tolerance, tolerance);
}

// test stepper getters for particle state
BOOST_AUTO_TEST_CASE(Getters) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magneticField->makeCache(magCtx),
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
      geoCtx, magneticField->makeCache(magCtx),
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  auto newPos4 = (pos4 + Vector4(1_mm, 2_mm, 3_mm, 20_ns)).eval();
  auto newPos = newPos4.segment<3>(ePos0);
  auto newTime = newPos4[eTime];
  auto newUnitDir = (unitDir + Vector3(1, -1, -1)).normalized();
  auto newAbsMom = 0.9 * absMom;

  // example surface and bound parameters at the updated position
  auto plane = Surface::makeShared<PlaneSurface>(newPos, newUnitDir);
  auto params = BoundTrackParameters::create(plane, geoCtx, newPos4, newUnitDir,
                                             charge, cov)
                    .value();
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
  BOOST_CHECK(params.covariance().has_value());
  stepper.update(state, freeParams, params.parameters(), *params.covariance(),
                 *plane);
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
      geoCtx, magneticField->makeCache(magCtx),
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  auto newPos = (pos + Vector3(1_mm, 2_mm, 3_mm)).eval();
  auto newTime = time + 20_ns;
  auto newUnitDir = (unitDir + Vector3(1, -1, -1)).normalized();
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
      geoCtx, magneticField->makeCache(magCtx),
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);
  // example surface at the current state position
  auto plane = Surface::makeShared<PlaneSurface>(pos, unitDir);

  auto&& [pars, jac, pathLength] = stepper.boundState(state, *plane).value();
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
      geoCtx, magneticField->makeCache(magCtx),
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
      geoCtx, magneticField->makeCache(magCtx),
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance));
  state.stepping.covTransport = false;

  // ensure step does not result in an error
  auto res = stepper.step(state, navigator);
  BOOST_CHECK(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_CHECK_EQUAL(state.stepping.stepSize.value(), navDir * stepSize);
  BOOST_CHECK_EQUAL(state.stepping.stepSize.value(), h);

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
      geoCtx, magneticField->makeCache(magCtx),
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance));
  state.stepping.covTransport = true;

  // ensure step does not result in an error
  auto res = stepper.step(state, navigator);
  BOOST_CHECK(res.ok());

  // extract the actual step size
  auto h = res.value();
  BOOST_CHECK_EQUAL(state.stepping.stepSize.value(), navDir * stepSize);
  BOOST_CHECK_EQUAL(state.stepping.stepSize.value(), h);

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

  stepper.transportCovarianceToCurvilinear(state.stepping);
  BOOST_CHECK_NE(state.stepping.cov, cov);
}

// test state reset method
BOOST_AUTO_TEST_CASE(Reset) {
  Stepper stepper(magneticField);
  MockPropagatorState state(Stepper::State(
      geoCtx, magneticField->makeCache(magCtx),
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance));
  state.stepping.covTransport = true;

  // ensure step does not result in an error
  stepper.step(state, navigator);

  // Construct the parameters
  Vector3 newPos(1.5, -2.5, 3.5);
  auto newAbsMom = 4.2 * absMom;
  double newTime = 7.5;
  double newCharge = 1.;
  BoundSymMatrix newCov = 8.5 * Covariance::Identity();
  CurvilinearTrackParameters cp(makeVector4(newPos, newTime), unitDir,
                                newAbsMom, newCharge, newCov);
  FreeVector freeParams = detail::transformBoundToFreeParameters(
      cp.referenceSurface(), geoCtx, cp.parameters());
  NavigationDirection ndir = NavigationDirection::Forward;
  double stepSize = -256.;

  auto copyState = [&](auto& field, const auto& other) {
    using field_t = std::decay_t<decltype(field)>;
    std::decay_t<decltype(other)> copy(geoCtx, field.makeCache(magCtx), cp,
                                       ndir, stepSize, tolerance);

    copy.state_ready = other.state_ready;
    copy.navDir = other.navDir;
    copy.useJacobian = other.useJacobian;
    copy.step = other.step;
    copy.maxPathLength = other.maxPathLength;
    copy.mcondition = other.mcondition;
    copy.needgradient = other.needgradient;
    copy.newfield = other.newfield;
    copy.field = other.field;
    copy.pVector = other.pVector;
    std::copy(std::begin(other.parameters), std::end(other.parameters),
              std::begin(copy.parameters));
    copy.covariance = other.covariance;
    copy.covTransport = other.covTransport;
    std::copy(std::begin(other.jacobian), std::end(other.jacobian),
              std::begin(copy.jacobian));
    copy.pathAccumulated = other.pathAccumulated;
    copy.stepSize = other.stepSize;
    copy.previousStepSize = other.previousStepSize;
    copy.tolerance = other.tolerance;

    copy.fieldCache =
        MagneticFieldProvider::Cache::make<typename field_t::Cache>(
            other.fieldCache.template get<typename field_t::Cache>());

    copy.geoContext = other.geoContext;
    copy.debug = other.debug;
    copy.debugString = other.debugString;
    copy.debugPfxWidth = other.debugPfxWidth;
    copy.debugMsgWidth = other.debugMsgWidth;

    return copy;
  };

  // Reset all possible parameters
  Stepper::State stateCopy(copyState(*magneticField, state.stepping));
  BOOST_CHECK(cp.covariance().has_value());
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface(), ndir, stepSize);
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, newCov);
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
  BOOST_CHECK_EQUAL(stateCopy.stepSize.value(), ndir * stepSize);
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize,
                    state.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(stateCopy.tolerance, state.stepping.tolerance);

  // Reset all possible parameters except the step size
  stateCopy = copyState(*magneticField, state.stepping);
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface(), ndir);
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, newCov);
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
  BOOST_CHECK_EQUAL(stateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize,
                    state.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(stateCopy.tolerance, state.stepping.tolerance);

  // Reset the least amount of parameters
  stateCopy = copyState(*magneticField, state.stepping);
  stepper.resetState(stateCopy, cp.parameters(), *cp.covariance(),
                     cp.referenceSurface());
  // Test all components
  BOOST_CHECK(stateCopy.covTransport);
  BOOST_CHECK_EQUAL(*stateCopy.covariance, newCov);
  BOOST_CHECK_EQUAL(stepper.position(stateCopy),
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(stepper.direction(stateCopy),
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(stepper.momentum(stateCopy),
                    std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(stepper.charge(stateCopy), stepper.charge(state.stepping));
  BOOST_CHECK_EQUAL(stepper.time(stateCopy), freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(stateCopy.navDir, NavigationDirection::Forward);
  BOOST_CHECK_EQUAL(stateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(stateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(stateCopy.previousStepSize,
                    state.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(stateCopy.tolerance, state.stepping.tolerance);

  // Reset using different surface shapes
  // 1) Disc surface
  // Setting some parameters
  newPos << 0.5, -1.5, 0.;
  newAbsMom *= 1.23;
  newTime = 8.4;
  newCharge = -1.;
  newCov = 10.9 * Covariance::Identity();
  Transform3 trafo = Transform3::Identity();
  auto disc = Surface::makeShared<DiscSurface>(trafo);
  auto boundDisc =
      BoundTrackParameters::create(disc, geoCtx, makeVector4(newPos, newTime),
                                   unitDir, newAbsMom, newCharge, newCov)
          .value();

  // Reset the state and test
  Stepper::State stateDisc = copyState(*magneticField, state.stepping);
  BOOST_CHECK(boundDisc.covariance().has_value());
  stepper.resetState(stateDisc, boundDisc.parameters(), *boundDisc.covariance(),
                     boundDisc.referenceSurface());

  CHECK_NE_COLLECTIONS(stateDisc.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(stateDisc.pVector, state.stepping.pVector);

  // 2) Perigee surface
  // Setting some parameters
  newPos << 1.5, -2.5, 3.5;
  newAbsMom *= 0.45;
  newTime = 2.3;
  newCharge = 1.;
  newCov = 8.7 * Covariance::Identity();
  auto perigee = Surface::makeShared<PerigeeSurface>(trafo);
  auto boundPerigee = BoundTrackParameters::create(
                          perigee, geoCtx, makeVector4(newPos, newTime),
                          unitDir, newAbsMom, newCharge, newCov)
                          .value();

  // Reset the state and test
  Stepper::State statePerigee = copyState(*magneticField, state.stepping);
  BOOST_CHECK(boundPerigee.covariance().has_value());
  stepper.resetState(statePerigee, boundPerigee.parameters(),
                     *boundPerigee.covariance(),
                     boundPerigee.referenceSurface());
  CHECK_NE_COLLECTIONS(statePerigee.pVector, stateCopy.pVector);
  CHECK_NE_COLLECTIONS(statePerigee.pVector, state.stepping.pVector);
  CHECK_NE_COLLECTIONS(statePerigee.pVector, stateDisc.pVector);

  // 3) Straw surface
  // Use the same parameters as for previous Perigee surface
  auto straw = Surface::makeShared<StrawSurface>(trafo);
  auto boundStraw =
      BoundTrackParameters::create(straw, geoCtx, makeVector4(newPos, newTime),
                                   unitDir, newAbsMom, newCharge, newCov)
          .value();

  // Reset the state and test
  Stepper::State stateStraw = copyState(*magneticField, state.stepping);
  BOOST_CHECK(boundStraw.covariance().has_value());
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
      geoCtx, magneticField->makeCache(magCtx),
      CurvilinearTrackParameters(pos4, unitDir, absMom, charge, cov), navDir,
      stepSize, tolerance);

  // TODO figure out why this fails and what it should be
  // BOOST_CHECK_EQUAL(stepper.overstepLimit(state), tolerance);

  stepper.setStepSize(state, 5_cm);
  BOOST_CHECK_EQUAL(state.previousStepSize, navDir * stepSize);
  BOOST_CHECK_EQUAL(state.stepSize.value(), 5_cm);

  stepper.releaseStepSize(state);
  BOOST_CHECK_EQUAL(state.stepSize.value(), navDir * stepSize);
}

// test step size modification with target surfaces
BOOST_AUTO_TEST_CASE(StepSizeSurface) {
  Stepper stepper(magneticField);
  Stepper::State state(
      geoCtx, magneticField->makeCache(magCtx),
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
  BOOST_CHECK_EQUAL(state.stepSize.value(), distance);

  // start with a different step size
  state.stepSize.setValue(navDir * stepSize);
  stepper.updateStepSize(
      state,
      target->intersect(state.geoContext, stepper.position(state),
                        state.navDir * stepper.direction(state), false),
      true);
  BOOST_CHECK_EQUAL(state.stepSize.value(), distance);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
