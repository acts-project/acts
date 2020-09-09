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
  PropState(StraightLineStepper::State sState) : stepping(sState) {}
  /// State of the straight line stepper
  StraightLineStepper::State stepping;
  /// Propagator options which only carry the particle's mass
  struct {
    double mass = 42.;
  } options;
};

static constexpr auto eps = 2 * std::numeric_limits<double>::epsilon();

/// These tests are aiming to test whether the state setup is working properly
BOOST_AUTO_TEST_CASE(straight_line_stepper_state_test) {
  // Set up some variables
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;

  Vector3D pos(1., 2., 3.);
  Vector3D dir(4., 5., 6.);
  double time = 7.;
  double absMom = 8.;
  double charge = -1.;

  // Test charged parameters without covariance matrix
  CurvilinearTrackParameters cp(makeVector4(pos, time), dir, absMom, charge);
  StraightLineStepper::State slsState(tgContext, mfContext, cp, ndir, stepSize,
                                      tolerance);

  // Test the result & compare with the input/test for reasonable members
  BOOST_CHECK_EQUAL(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());
  BOOST_CHECK(!slsState.covTransport);
  BOOST_CHECK_EQUAL(slsState.cov, Covariance::Zero());
  CHECK_CLOSE_OR_SMALL(slsState.pos, pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(slsState.dir, dir.normalized(), eps, eps);
  CHECK_CLOSE_REL(slsState.p, absMom, eps);
  BOOST_CHECK_EQUAL(slsState.q, charge);
  CHECK_CLOSE_OR_SMALL(slsState.t, time, eps, eps);
  BOOST_CHECK_EQUAL(slsState.navDir, ndir);
  BOOST_CHECK_EQUAL(slsState.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsState.stepSize, ndir * stepSize);
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
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;

  // Construct the parameters
  Vector3D pos(1., 2., 3.);
  Vector3D dir = Vector3D(4., 5., 6.).normalized();
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
  BOOST_CHECK_EQUAL(sls.position(slsState), slsState.pos);
  BOOST_CHECK_EQUAL(sls.direction(slsState), slsState.dir);
  BOOST_CHECK_EQUAL(sls.momentum(slsState), slsState.p);
  BOOST_CHECK_EQUAL(sls.charge(slsState), slsState.q);
  BOOST_CHECK_EQUAL(sls.time(slsState), slsState.t);

  //~ BOOST_CHECK_EQUAL(sls.overstepLimit(slsState), tolerance);

  // Step size modifies
  const std::string originalStepSize = slsState.stepSize.toString();

  sls.setStepSize(slsState, 1337.);
  BOOST_CHECK_EQUAL(slsState.previousStepSize, ndir * stepSize);
  BOOST_CHECK_EQUAL(slsState.stepSize, 1337.);

  sls.releaseStepSize(slsState);
  BOOST_CHECK_EQUAL(slsState.stepSize, -123.);
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
  Vector3D newPos(2., 4., 8.);
  Vector3D newMom(3., 9., 27.);
  double newTime(321.);
  sls.update(slsState, newPos, newMom.normalized(), newMom.norm(), newTime);
  BOOST_CHECK_EQUAL(slsState.pos, newPos);
  BOOST_CHECK_EQUAL(slsState.dir, newMom.normalized());
  BOOST_CHECK_EQUAL(slsState.p, newMom.norm());
  BOOST_CHECK_EQUAL(slsState.q, charge);
  BOOST_CHECK_EQUAL(slsState.t, newTime);

  // The covariance transport
  slsState.cov = cov;
  sls.covarianceTransport(slsState);
  BOOST_CHECK_NE(slsState.cov, cov);
  BOOST_CHECK_NE(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());

  // Perform a step without and with covariance transport
  slsState.cov = cov;
  PropState ps(slsState);

  ps.stepping.covTransport = false;
  double h = sls.step(ps).value();
  BOOST_CHECK_EQUAL(ps.stepping.stepSize, ndir * stepSize);
  BOOST_CHECK_EQUAL(ps.stepping.stepSize, h);
  CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
  BOOST_CHECK_GT(ps.stepping.pos.norm(), newPos.norm());
  BOOST_CHECK_EQUAL(ps.stepping.dir, newMom.normalized());
  BOOST_CHECK_EQUAL(ps.stepping.p, newMom.norm());
  BOOST_CHECK_EQUAL(ps.stepping.q, charge);
  BOOST_CHECK_LT(ps.stepping.t, newTime);
  BOOST_CHECK_EQUAL(ps.stepping.derivative, FreeVector::Zero());
  BOOST_CHECK_EQUAL(ps.stepping.jacTransport, FreeMatrix::Identity());

  ps.stepping.covTransport = true;
  double h2 = sls.step(ps).value();
  BOOST_CHECK_EQUAL(ps.stepping.stepSize, ndir * stepSize);
  BOOST_CHECK_EQUAL(h2, h);
  CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
  BOOST_CHECK_GT(ps.stepping.pos.norm(), newPos.norm());
  BOOST_CHECK_EQUAL(ps.stepping.dir, newMom.normalized());
  BOOST_CHECK_EQUAL(ps.stepping.p, newMom.norm());
  BOOST_CHECK_EQUAL(ps.stepping.q, charge);
  BOOST_CHECK_LT(ps.stepping.t, newTime);
  BOOST_CHECK_NE(ps.stepping.derivative, FreeVector::Zero());
  BOOST_CHECK_NE(ps.stepping.jacTransport, FreeMatrix::Identity());

  /// Test the state reset
  // Construct the parameters
  Vector3D pos2(1.5, -2.5, 3.5);
  Vector3D dir2 = Vector3D(4.5, -5.5, 6.5).normalized();
  double time2 = 7.5;
  double absMom2 = 8.5;
  double charge2 = 1.;
  BoundSymMatrix cov2 = 8.5 * Covariance::Identity();
  CurvilinearTrackParameters cp2(makeVector4(pos2, time2), dir2, absMom2,
                                 charge2, cov2);
  FreeVector freeParams = detail::transformBoundToFreeParameters(
      cp2.referenceSurface(), tgContext, cp2.parameters());
  ndir = forward;
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
  BOOST_CHECK_EQUAL(slsStateCopy.pos,
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(slsStateCopy.dir,
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(slsStateCopy.p, std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(slsStateCopy.q, ps.stepping.q);
  BOOST_CHECK_EQUAL(slsStateCopy.t, freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(slsStateCopy.navDir, ndir);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize, ndir * stepSize2);
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
  BOOST_CHECK_EQUAL(slsStateCopy.pos,
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(slsStateCopy.dir,
                    freeParams.template segment<3>(eFreeDir0));
  BOOST_CHECK_EQUAL(slsStateCopy.p, std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(slsStateCopy.q, ps.stepping.q);
  BOOST_CHECK_EQUAL(slsStateCopy.t, freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(slsStateCopy.navDir, ndir);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize,
                    ndir * std::numeric_limits<double>::max());
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
  BOOST_CHECK_EQUAL(slsStateCopy.pos,
                    freeParams.template segment<3>(eFreePos0));
  BOOST_CHECK_EQUAL(slsStateCopy.dir,
                    freeParams.template segment<3>(eFreeDir0).normalized());
  BOOST_CHECK_EQUAL(slsStateCopy.p, std::abs(1. / freeParams[eFreeQOverP]));
  BOOST_CHECK_EQUAL(slsStateCopy.q, ps.stepping.q);
  BOOST_CHECK_EQUAL(slsStateCopy.t, freeParams[eFreeTime]);
  BOOST_CHECK_EQUAL(slsStateCopy.navDir, forward);
  BOOST_CHECK_EQUAL(slsStateCopy.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(slsStateCopy.stepSize, std::numeric_limits<double>::max());
  BOOST_CHECK_EQUAL(slsStateCopy.previousStepSize,
                    ps.stepping.previousStepSize);
  BOOST_CHECK_EQUAL(slsStateCopy.tolerance, ps.stepping.tolerance);

  /// Repeat with surface related methods
  auto plane = Surface::makeShared<PlaneSurface>(pos, dir);
  BoundTrackParameters bp(plane, tgContext, makeVector4(pos, time), dir,
                          charge / absMom, cov);
  slsState = StraightLineStepper::State(tgContext, mfContext, cp, ndir,
                                        stepSize, tolerance);

  // Test the intersection in the context of a surface
  auto targetSurface =
      Surface::makeShared<PlaneSurface>(pos + ndir * 2. * dir, dir);
  sls.updateSurfaceStatus(slsState, *targetSurface, BoundaryCheck(false));
  CHECK_CLOSE_ABS(slsState.stepSize.value(ConstrainedStep::actor), ndir * 2.,
                  1e-6);

  // Test the step size modification in the context of a surface
  sls.updateStepSize(
      slsState,
      targetSurface->intersect(slsState.geoContext, slsState.pos,
                               slsState.navDir * slsState.dir, false),
      false);
  CHECK_CLOSE_ABS(slsState.stepSize, 2, 1e-6);
  slsState.stepSize = ndir * stepSize;
  sls.updateStepSize(
      slsState,
      targetSurface->intersect(slsState.geoContext, slsState.pos,
                               slsState.navDir * slsState.dir, false),
      true);
  CHECK_CLOSE_ABS(slsState.stepSize, 2, 1e-6);

  // Test the bound state construction
  auto boundState = sls.boundState(slsState, *plane);
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
  sls.covarianceTransport(slsState, *plane);
  BOOST_CHECK_NE(slsState.cov, cov);
  BOOST_CHECK_NE(slsState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(slsState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(slsState.derivative, FreeVector::Zero());

  // Update in context of a surface
  freeParams = detail::transformBoundToFreeParameters(
      bp.referenceSurface(), tgContext, bp.parameters());
  freeParams.segment<3>(eFreePos0) *= 2;
  freeParams[eFreeTime] *= 2;
  freeParams.segment<3>(eFreeDir0) *= 2;
  freeParams[eFreeQOverP] *= -0.5;

  sls.update(slsState, freeParams, 2 * (*bp.covariance()));
  CHECK_CLOSE_OR_SMALL(slsState.pos, 2. * pos, eps, eps);
  CHECK_CLOSE_OR_SMALL(slsState.dir, dir, eps, eps);
  CHECK_CLOSE_REL(slsState.p, 2. * absMom, eps);
  BOOST_CHECK_EQUAL(slsState.q, 1. * charge);
  CHECK_CLOSE_OR_SMALL(slsState.t, 2. * time, eps, eps);
  CHECK_CLOSE_COVARIANCE(slsState.cov, Covariance(2. * cov), 1e-6);
}
}  // namespace Test
}  // namespace Acts
