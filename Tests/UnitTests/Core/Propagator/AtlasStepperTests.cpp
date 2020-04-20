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

namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Jacobian = BoundMatrix;

/// @brief Simplified propagator state
struct PropState {
  /// @brief Constructor
  PropState(AtlasStepper<ConstantBField>::State sState) : stepping(sState) {}
  /// State of the Atlas stepper
  AtlasStepper<ConstantBField>::State stepping;
  /// Propagator options which only carry the relevant components
  struct {
    double mass = 42.;
    double tolerance = 1e-4;
  } options;
};

/// These tests are aiming to test whether the state setup is working properly
BOOST_AUTO_TEST_CASE(atlas_stepper_state_test) {
  // Set up some variables
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;
  ConstantBField bField(Vector3D(1., 2.5, 33.33));

  Vector3D pos(1., 2., 3.);
  Vector3D mom(4., 5., 6.);
  Vector3D dir = mom.normalized();
  double time = 7.;
  double charge = -1.;

  // Test charged parameters without covariance matrix
  CurvilinearParameters cp(std::nullopt, pos, mom, charge, time);
  AtlasStepper<ConstantBField>::State asState(tgContext, mfContext, cp, ndir,
                                              stepSize, tolerance);

  // Test the result & compare with the input/test for reasonable members
  BOOST_TEST(!asState.covTransport);
  BOOST_TEST(asState.covariance == nullptr);
  BOOST_TEST(asState.pVector[0] == pos.x());
  BOOST_TEST(asState.pVector[1] == pos.y());
  BOOST_TEST(asState.pVector[2] == pos.z());
  BOOST_TEST(asState.pVector[3] == time);
  CHECK_CLOSE_ABS(asState.pVector[4], dir.x(), 1e-6);
  CHECK_CLOSE_ABS(asState.pVector[5], dir.y(), 1e-6);
  CHECK_CLOSE_ABS(asState.pVector[6], dir.z(), 1e-6);
  BOOST_TEST(asState.pVector[7] == charge / mom.norm());
  BOOST_TEST(asState.navDir == ndir);
  BOOST_TEST(asState.pathAccumulated == 0.);
  BOOST_TEST(asState.stepSize == ndir * stepSize);
  BOOST_TEST(asState.previousStepSize == 0.);
  BOOST_TEST(asState.tolerance == tolerance);

  // Test with covariance matrix
  Covariance cov = 8. * Covariance::Identity();
  cp = CurvilinearParameters(cov, pos, mom, charge, time);
  asState = AtlasStepper<ConstantBField>::State(tgContext, mfContext, cp, ndir,
                                                stepSize, tolerance);
  BOOST_TEST(asState.covTransport);
  BOOST_TEST(*asState.covariance == cov);
}

/// These tests are aiming to test the functions of the AtlasStepper>
/// The numerical correctness of the stepper is tested in the integration tests
BOOST_AUTO_TEST_CASE(atlas_stepper_test) {
  // Set up some variables for the state
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;
  ConstantBField bField(Vector3D(1., 2.5, 33.33));

  // Construct the parameters
  Vector3D pos(1., 2., 3.);
  Vector3D mom(4., 5., 6.);
  double time = 7.;
  double charge = -1.;
  Covariance cov = 8. * Covariance::Identity();
  CurvilinearParameters cp(cov, pos, mom, charge, time);

  // Build the state and the stepper
  AtlasStepper<ConstantBField>::State asState(tgContext, mfContext, cp, ndir,
                                              stepSize, tolerance);
  AtlasStepper<ConstantBField> as;

  // Test the getters
  BOOST_TEST(as.position(asState) == pos);
  CHECK_CLOSE_ABS(as.direction(asState), mom.normalized(), 1e-6);
  BOOST_TEST(as.momentum(asState) == mom.norm());
  BOOST_TEST(as.charge(asState) == charge);
  BOOST_TEST(as.time(asState) == time);

  // BOOST_TEST(as.overstepLimit(asState) == tolerance);

  // Step size modifies
  const std::string originalStepSize = asState.stepSize.toString();

  as.setStepSize(asState, 1337.);
  BOOST_TEST(asState.previousStepSize == ndir * stepSize);
  BOOST_TEST(asState.stepSize == 1337.);

  as.releaseStepSize(asState);
  BOOST_TEST(asState.stepSize == -123.);
  BOOST_TEST(as.outputStepSize(asState) == originalStepSize);

  // Test the curvilinear state construction
  auto curvState = as.curvilinearState(asState);
  auto curvPars = std::get<0>(curvState);
  CHECK_CLOSE_ABS(curvPars.position(), cp.position(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.momentum(), cp.momentum(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.charge(), cp.charge(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.time(), cp.time(), 1e-6);
  BOOST_TEST(curvPars.covariance().has_value());
  BOOST_TEST(*curvPars.covariance() != cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(curvState),
                         BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_ABS(std::get<2>(curvState), 0., 1e-6);

  // Test the update method
  Vector3D newPos(2., 4., 8.);
  Vector3D newMom(3., 9., 27.);
  double newTime(321.);
  as.update(asState, newPos, newMom.normalized(), newMom.norm(), newTime);
  BOOST_TEST(as.position(asState) == newPos);
  BOOST_TEST(as.direction(asState) == newMom.normalized());
  BOOST_TEST(as.momentum(asState) == newMom.norm());
  BOOST_TEST(as.charge(asState) == charge);
  BOOST_TEST(as.time(asState) == newTime);

  // The covariance transport
  asState.cov = cov;
  as.covarianceTransport(asState);
  BOOST_TEST(asState.cov != cov);

  // Perform a step without and with covariance transport
  asState.cov = cov;
  PropState ps(asState);

  ps.stepping.covTransport = false;
  double h = as.step(ps).value();
  BOOST_TEST(ps.stepping.stepSize == ndir * stepSize);
  BOOST_TEST(ps.stepping.stepSize == h);
  CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
  BOOST_TEST(as.position(ps.stepping).norm() > newPos.norm());
  CHECK_CLOSE_ABS(as.direction(ps.stepping), newMom.normalized(), 1e-6);
  BOOST_TEST(as.momentum(ps.stepping) == newMom.norm());
  BOOST_TEST(as.charge(ps.stepping) == charge);
  BOOST_TEST(as.time(ps.stepping) < newTime);

  ps.stepping.covTransport = true;
  double h2 = as.step(ps).value();
  BOOST_TEST(ps.stepping.stepSize == ndir * stepSize);
  BOOST_TEST(h2 == h);
  CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
  BOOST_TEST(as.position(ps.stepping).norm() > newPos.norm());
  CHECK_CLOSE_ABS(as.direction(ps.stepping), newMom.normalized(), 1e-6);
  BOOST_TEST(as.momentum(ps.stepping) == newMom.norm());
  BOOST_TEST(as.charge(ps.stepping) == charge);
  BOOST_TEST(as.time(ps.stepping) < newTime);

  /// Repeat with surface related methods
  auto plane = Surface::makeShared<PlaneSurface>(pos, mom.normalized());
  BoundParameters bp(tgContext, cov, pos, mom, charge, time, plane);
  asState = AtlasStepper<ConstantBField>::State(tgContext, mfContext, cp, ndir,
                                                stepSize, tolerance);

  // Test the intersection in the context of a surface
  auto targetSurface = Surface::makeShared<PlaneSurface>(
      pos + ndir * 2. * mom.normalized(), mom.normalized());
  as.updateSurfaceStatus(asState, *targetSurface, BoundaryCheck(false));
  BOOST_TEST(asState.stepSize.value(ConstrainedStep::actor), ndir * 2.);

  // Test the step size modification in the context of a surface
  as.updateStepSize(
      asState,
      targetSurface->intersect(asState.geoContext, as.position(asState),
                               asState.navDir * as.direction(asState), false),
      false);
  BOOST_TEST(asState.stepSize == 2.);
  asState.stepSize = ndir * stepSize;
  as.updateStepSize(
      asState,
      targetSurface->intersect(asState.geoContext, as.position(asState),
                               asState.navDir * as.direction(asState), false),
      true);
  BOOST_TEST(asState.stepSize == 2.);

  // Test the bound state construction
  auto boundState = as.boundState(asState, *plane);
  auto boundPars = std::get<0>(boundState);
  CHECK_CLOSE_ABS(boundPars.position(), bp.position(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.momentum(), bp.momentum(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.charge(), bp.charge(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.time(), bp.time(), 1e-6);
  BOOST_TEST(boundPars.covariance().has_value());
  BOOST_TEST(*boundPars.covariance() != cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(boundState),
                         BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_ABS(std::get<2>(boundState), 0., 1e-6);

  // Update in context of a surface
  BoundParameters bpTarget(tgContext, 2. * cov, 2. * pos, 2. * mom,
                           -1. * charge, 2. * time, targetSurface);
  as.update(asState, bpTarget);
  BOOST_TEST(as.position(asState) == 2. * pos);
  CHECK_CLOSE_ABS(as.direction(asState), mom.normalized(), 1e-6);
  BOOST_TEST(as.momentum(asState) == 2. * mom.norm());
  BOOST_TEST(as.charge(asState) == -1. * charge);
  BOOST_TEST(as.time(asState) == 2. * time);
  CHECK_CLOSE_COVARIANCE(*asState.covariance, Covariance(2. * cov), 1e-6);

  // Transport the covariance in the context of a surface
  as.covarianceTransport(asState, *plane);
  BOOST_TEST(asState.cov != cov);
}
}  // namespace Test
}  // namespace Acts