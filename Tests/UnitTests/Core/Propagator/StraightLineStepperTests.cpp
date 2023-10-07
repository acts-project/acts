// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <limits>

namespace tt = boost::test_tools;
using Acts::VectorHelpers::makeVector4;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Jacobian = BoundMatrix;

/// @brief Simplified propagator state
struct PropState {
  /// @brief Constructor
  PropState(StraightLineStepper::State sState) : stepping(std::move(sState)) {}
  /// State of the straight line stepper
  StraightLineStepper::State stepping;
  /// Propagator options which only carry the particle's mass
  struct {
    double mass = 42.;
  } options;
};

struct MockNavigator {};

static constexpr MockNavigator mockNavigator;

static constexpr auto eps = 2 * std::numeric_limits<double>::epsilon();

/// These tests are aiming to test whether the state setup is working properly
BOOST_AUTO_TEST_CASE(straight_line_stepper_state_test) {
  // Set up some variables
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  NavigationDirection ndir = NavigationDirection::Backward;
  double stepSize = 123.;
  double tolerance = 234.;

  Vector3 pos(1., 2., 3.);
  Vector3 dir(4., 5., 6.);
  double time = 7.;
  double absMom = 8.;
  double charge = -1.;

  // Test charged parameters without covariance matrix
  CurvilinearTrackParameters cp(makeVector4(pos, time), dir, absMom, charge);
  StraightLineStepper::State slsState(tgContext, mfContext, cp, ndir, stepSize,
                                      tolerance);

  StraightLineStepper sls;

  // Test the result & compare with the input/test for reasonable members
  BOOST_CHECK_EQUAL(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());
  BOOST_CHECK(!slsState.covTransport);
  BOOST_CHECK_EQUAL(slsState.cov, Covariance::Zero());
  CHECK_CLOSE_OR_SMALL(sls.position(slsState), pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(sls.direction(slsState), dir.normalized(), eps, eps);
  CHECK_CLOSE_REL(sls.momentum(slsState), absMom, eps);
  BOOST_CHECK_EQUAL(sls.charge(slsState), charge);
  CHECK_CLOSE_OR_SMALL(sls.time(slsState), time, eps, eps);
  BOOST_CHECK_EQUAL(slsState.navDir, ndir);
  BOOST_CHECK_EQUAL(slsState.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), ndir * stepSize);
  BOOST_CHECK_EQUAL(slsState.previousStepSize, 0.);
  BOOST_CHECK_EQUAL(slsState.tolerance, tolerance);

  // Test without charge and covariance matrix
  NeutralCurvilinearTrackParameters ncp(makeVector4(pos, time), dir,
                                        1 / absMom);
  slsState = StraightLineStepper::State(tgContext, mfContext, ncp, ndir,
                                        stepSize, tolerance);
  BOOST_CHECK_EQUAL(slsState.q, 0.);

  // Test with covariance matrix
  Covariance cov = 8. * Covariance::Identity();
  ncp = NeutralCurvilinearTrackParameters(makeVector4(pos, time), dir,
                                          1 / absMom, cov);
  slsState = StraightLineStepper::State(tgContext, mfContext, ncp, ndir,
                                        stepSize, tolerance);
  BOOST_CHECK_NE(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK(slsState.covTransport);
  BOOST_CHECK_EQUAL(slsState.cov, cov);
}

/// These tests are aiming to test the functions of the StraightLineStepper
/// The numerical correctness of the stepper is tested in the integration tests
BOOST_AUTO_TEST_CASE(straight_line_stepper_test) {
  // Set up some variables for the state
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  NavigationDirection ndir = NavigationDirection::Backward;
  double stepSize = 123.;
  double tolerance = 234.;

  // Construct the parameters
  Vector3 pos(1., 2., 3.);
  Vector3 dir = Vector3(4., 5., 6.).normalized();
  double time = 7.;
  double absMom = 8.;
  double charge = -1.;
  Covariance cov = 8. * Covariance::Identity();
  CurvilinearTrackParameters cp(makeVector4(pos, time), dir, charge / absMom,
                                cov);

  // Build the state and the stepper
  StraightLineStepper::State slsState(tgContext, mfContext, cp, ndir, stepSize,
                                      tolerance);
  StraightLineStepper sls;

  // Test the getters
  CHECK_CLOSE_ABS(sls.position(slsState), pos, 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsState), dir, 1e-6);
  CHECK_CLOSE_ABS(sls.momentum(slsState), absMom, 1e-6);
  BOOST_CHECK_EQUAL(sls.charge(slsState), charge);
  BOOST_CHECK_EQUAL(sls.time(slsState), time);

  //~ BOOST_CHECK_EQUAL(sls.overstepLimit(slsState), tolerance);

  // Step size modifies
  const std::string originalStepSize = slsState.stepSize.toString();

  sls.setStepSize(slsState, 1337.);
  BOOST_CHECK_EQUAL(slsState.previousStepSize, ndir * stepSize);
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), 1337.);

  sls.releaseStepSize(slsState);
  BOOST_CHECK_EQUAL(slsState.stepSize.value(), -123.);
  BOOST_CHECK_EQUAL(sls.outputStepSize(slsState), originalStepSize);

  // Test the curvilinear state construction
  auto curvState = sls.curvilinearState(slsState);
  auto curvPars = std::get<0>(curvState);
  CHECK_CLOSE_ABS(curvPars.position(tgContext), cp.position(tgContext), 1e-6);
  CHECK_CLOSE_ABS(curvPars.momentum(), cp.momentum(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.charge(), cp.charge(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.time(), cp.time(), 1e-6);
  BOOST_CHECK(curvPars.covariance().has_value());
  BOOST_CHECK_NE(*curvPars.covariance(), cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(curvState),
                         BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_ABS(std::get<2>(curvState), 0., 1e-6);

  // Test the update method
  Vector3 newPos(2., 4., 8.);
  Vector3 newMom(3., 9., 27.);
  double newTime(321.);
  sls.update(slsState, newPos, newMom.normalized(), newMom.norm(), newTime);
  CHECK_CLOSE_ABS(sls.position(slsState), newPos, 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsState), newMom.normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.momentum(slsState), newMom.norm(), 1e-6);
  BOOST_CHECK_EQUAL(sls.charge(slsState), charge);
  BOOST_CHECK_EQUAL(sls.time(slsState), newTime);

  // The covariance transport
  slsState.cov = cov;
  sls.transportCovarianceToCurvilinear(slsState);
  BOOST_CHECK_NE(slsState.cov, cov);
  BOOST_CHECK_NE(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());

  // Perform a step without and with covariance transport
  slsState.cov = cov;
  PropState ps(slsState);

  ps.stepping.covTransport = false;
  double h = sls.step(ps, mockNavigator).value();
  BOOST_CHECK_EQUAL(ps.stepping.stepSize.value(), ndir * stepSize);
  BOOST_CHECK_EQUAL(ps.stepping.stepSize.value(), h);
  CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
  BOOST_CHECK_GT(sls.position(ps.stepping).norm(), newPos.norm());
  CHECK_CLOSE_ABS(sls.direction(ps.stepping), newMom.normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.momentum(ps.stepping), newMom.norm(), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(ps.stepping), charge, 1e-6);
  BOOST_CHECK_LT(sls.time(ps.stepping), newTime);
  BOOST_CHECK_EQUAL(ps.stepping.derivative, FreeVector::Zero());
  BOOST_CHECK_EQUAL(ps.stepping.jacTransport, FreeMatrix::Identity());

  ps.stepping.covTransport = true;
  double h2 = sls.step(ps, mockNavigator).value();
  BOOST_CHECK_EQUAL(ps.stepping.stepSize.value(), ndir * stepSize);
  BOOST_CHECK_EQUAL(h2, h);
  CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
  BOOST_CHECK_GT(sls.position(ps.stepping).norm(), newPos.norm());
  CHECK_CLOSE_ABS(sls.direction(ps.stepping), newMom.normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.momentum(ps.stepping), newMom.norm(), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(ps.stepping), charge, 1e-6);
  BOOST_CHECK_LT(sls.time(ps.stepping), newTime);
  BOOST_CHECK_NE(ps.stepping.derivative, FreeVector::Zero());
  BOOST_CHECK_NE(ps.stepping.jacTransport, FreeMatrix::Identity());

  /// Test the state reset
  // Construct the parameters
  Vector3 pos2(1.5, -2.5, 3.5);
  Vector3 dir2 = Vector3(4.5, -5.5, 6.5).normalized();
  double time2 = 7.5;
  double absMom2 = 8.5;
  double charge2 = 1.;
  BoundSymMatrix cov2 = 8.5 * Covariance::Identity();
  CurvilinearTrackParameters cp2(makeVector4(pos2, time2), dir2, absMom2,
                                 charge2, cov2);
  BOOST_CHECK(cp2.covariance().has_value());
  FreeVector freeParams = detail::transformBoundToFreeParameters(
      cp2.referenceSurface(), tgContext, cp2.parameters());
  ndir = NavigationDirection::Forward;
  double stepSize2 = -2. * stepSize;

  // Reset all possible parameters
  StraightLineStepper::State slsStateCopy(ps.stepping);
  sls.resetState(slsStateCopy, cp2.parameters(), *cp2.covariance(),
                 cp2.referenceSurface(), ndir, stepSize2);
  // Test all components
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, ps.stepping.jacToGlobal);
  BOOST_CHECK_EQUAL(slsStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(slsStateCopy.covTransport);
  BOOST_CHECK_EQUAL(slsStateCopy.cov, cov2);
  CHECK_CLOSE_ABS(sls.position(slsStateCopy),
                  freeParams.template segment<3>(eFreePos0), 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsStateCopy),
                  freeParams.template segment<3>(eFreeDir0).normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.momentum(slsStateCopy),
                  std::abs(1. / freeParams[eFreeQOverP]), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(slsStateCopy), sls.charge(ps.stepping), 1e-6);
  CHECK_CLOSE_ABS(sls.time(slsStateCopy), freeParams[eFreeTime], 1e-6);
  BOOST_CHECK_EQUAL(slsStateCopy.navDir, ndir);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize.value(), ndir * stepSize2);
  BOOST_CHECK_EQUAL(slsStateCopy.previousStepSize,
                    ps.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(slsStateCopy.tolerance, ps.stepping.tolerance);

  // Reset all possible parameters except the step size
  slsStateCopy = ps.stepping;
  sls.resetState(slsStateCopy, cp2.parameters(), *cp2.covariance(),
                 cp2.referenceSurface(), ndir);
  // Test all components
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, ps.stepping.jacToGlobal);
  BOOST_CHECK_EQUAL(slsStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(slsStateCopy.covTransport);
  BOOST_CHECK_EQUAL(slsStateCopy.cov, cov2);
  CHECK_CLOSE_ABS(sls.position(slsStateCopy),
                  freeParams.template segment<3>(eFreePos0), 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsStateCopy),
                  freeParams.template segment<3>(eFreeDir0), 1e-6);
  CHECK_CLOSE_ABS(sls.momentum(slsStateCopy),
                  std::abs(1. / freeParams[eFreeQOverP]), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(slsStateCopy), sls.charge(ps.stepping), 1e-6);
  CHECK_CLOSE_ABS(sls.time(slsStateCopy), freeParams[eFreeTime], 1e-6);
  BOOST_CHECK_EQUAL(slsStateCopy.navDir, ndir);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(slsStateCopy.previousStepSize,
                    ps.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(slsStateCopy.tolerance, ps.stepping.tolerance);

  // Reset the least amount of parameters
  slsStateCopy = ps.stepping;
  sls.resetState(slsStateCopy, cp2.parameters(), *cp2.covariance(),
                 cp2.referenceSurface());
  // Test all components
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_NE(slsStateCopy.jacToGlobal, ps.stepping.jacToGlobal);
  BOOST_CHECK_EQUAL(slsStateCopy.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsStateCopy.derivative, FreeVector::Zero());
  BOOST_CHECK(slsStateCopy.covTransport);
  BOOST_CHECK_EQUAL(slsStateCopy.cov, cov2);
  CHECK_CLOSE_ABS(sls.position(slsStateCopy),
                  freeParams.template segment<3>(eFreePos0), 1e-6);
  CHECK_CLOSE_ABS(sls.direction(slsStateCopy),
                  freeParams.template segment<3>(eFreeDir0).normalized(), 1e-6);
  CHECK_CLOSE_ABS(sls.momentum(slsStateCopy),
                  std::abs(1. / freeParams[eFreeQOverP]), 1e-6);
  CHECK_CLOSE_ABS(sls.charge(slsStateCopy), sls.charge(ps.stepping), 1e-6);
  CHECK_CLOSE_ABS(sls.time(slsStateCopy), freeParams[eFreeTime], 1e-6);
  BOOST_CHECK_EQUAL(slsStateCopy.navDir, NavigationDirection::Forward);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize.value(),
                    std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(slsStateCopy.previousStepSize,
                    ps.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(slsStateCopy.tolerance, ps.stepping.tolerance);

  /// Repeat with surface related methods
  auto plane = Surface::makeShared<PlaneSurface>(pos, dir);
  auto bp =
      BoundTrackParameters::create(plane, tgContext, makeVector4(pos, time),
                                   dir, charge / absMom, cov)
          .value();
  slsState = StraightLineStepper::State(tgContext, mfContext, cp, ndir,
                                        stepSize, tolerance);

  // Test the intersection in the context of a surface
  auto targetSurface =
      Surface::makeShared<PlaneSurface>(pos + ndir * 2. * dir, dir);
  sls.updateSurfaceStatus(slsState, *targetSurface, BoundaryCheck(false));
  CHECK_CLOSE_ABS(slsState.stepSize.value(ConstrainedStep::actor), ndir * 2.,
                  1e-6);

  // Test the step size modification in the context of a surface
  sls.updateStepSize(slsState,
                     targetSurface->intersect(
                         slsState.geoContext, sls.position(slsState),
                         slsState.navDir * sls.direction(slsState), false),
                     false);
  CHECK_CLOSE_ABS(slsState.stepSize.value(), 2, 1e-6);
  slsState.stepSize.setValue(ndir * stepSize);
  sls.updateStepSize(slsState,
                     targetSurface->intersect(
                         slsState.geoContext, sls.position(slsState),
                         slsState.navDir * sls.direction(slsState), false),
                     true);
  CHECK_CLOSE_ABS(slsState.stepSize.value(), 2, 1e-6);

  // Test the bound state construction
  FreeToBoundCorrection freeToBoundCorrection(false);
  auto boundState =
      sls.boundState(slsState, *plane, true, freeToBoundCorrection).value();
  auto boundPars = std::get<0>(boundState);
  CHECK_CLOSE_ABS(boundPars.position(tgContext), bp.position(tgContext), 1e-6);
  CHECK_CLOSE_ABS(boundPars.momentum(), bp.momentum(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.charge(), bp.charge(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.time(), bp.time(), 1e-6);
  BOOST_CHECK(boundPars.covariance().has_value());
  BOOST_CHECK_NE(*boundPars.covariance(), cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(boundState),
                         BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_ABS(std::get<2>(boundState), 0., 1e-6);

  // Transport the covariance in the context of a surface
  sls.transportCovarianceToBound(slsState, *plane, freeToBoundCorrection);
  BOOST_CHECK_NE(slsState.cov, cov);
  BOOST_CHECK_NE(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());

  // Update in context of a surface
  freeParams = detail::transformBoundToFreeParameters(
      bp.referenceSurface(), tgContext, bp.parameters());
  freeParams.segment<3>(eFreePos0) *= 2;
  freeParams[eFreeTime] *= 2;

  BOOST_CHECK(bp.covariance().has_value());
  sls.update(slsState, freeParams, bp.parameters(), 2 * (*bp.covariance()),
             *plane);
  CHECK_CLOSE_OR_SMALL(sls.position(slsState), 2. * pos, eps, eps);
  BOOST_CHECK_EQUAL(sls.charge(slsState), 1. * charge);
  CHECK_CLOSE_OR_SMALL(sls.time(slsState), 2. * time, eps, eps);
  CHECK_CLOSE_COVARIANCE(slsState.cov, Covariance(2. * cov), 1e-6);
}
}  // namespace Test
}  // namespace Acts
