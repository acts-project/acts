// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Propagator/MultiEigenStepperLoop.hpp"

BOOST_AUTO_TEST_CASE(eigen_stepper_state_test) {
  // Set up some variables
/*
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;
  auto bField = std::make_shared<ConstantBField>(Vector3(1., 2.5, 33.33));

  Vector3 pos(1., 2., 3.);
  Vector3 dir(4., 5., 6.);
  double time = 7.;
  double absMom = 8.;
  double charge = -1.;

  // Test charged parameters without covariance matrix
  CurvilinearTrackParameters cp(makeVector4(pos, time), dir, charge / absMom);
  EigenStepper<>::State esState(tgContext, bField->makeCache(mfContext), cp,
                                ndir, stepSize, tolerance);

  EigenStepper<> es(bField);

  // Test the result & compare with the input/test for reasonable members
  BOOST_CHECK_EQUAL(esState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK_EQUAL(esState.jacTransport, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(esState.derivative, FreeVector::Zero());
  BOOST_CHECK(!esState.covTransport);
  BOOST_CHECK_EQUAL(esState.cov, Covariance::Zero());
  BOOST_CHECK_EQUAL(esState.navDir, ndir);
  BOOST_CHECK_EQUAL(esState.pathAccumulated, 0.);
  BOOST_CHECK_EQUAL(esState.stepSize, ndir * stepSize);
  BOOST_CHECK_EQUAL(esState.previousStepSize, 0.);
  BOOST_CHECK_EQUAL(esState.tolerance, tolerance);

  // Test without charge and covariance matrix
  NeutralCurvilinearTrackParameters ncp(makeVector4(pos, time), dir,
                                        1 / absMom);
  esState = EigenStepper<>::State(tgContext, bField->makeCache(mfContext), ncp,
                                  ndir, stepSize, tolerance);
  BOOST_CHECK_EQUAL(es.charge(esState), 0.);

  // Test with covariance matrix
  Covariance cov = 8. * Covariance::Identity();
  ncp = NeutralCurvilinearTrackParameters(makeVector4(pos, time), dir,
                                          1 / absMom, cov);
  esState = EigenStepper<>::State(tgContext, bField->makeCache(mfContext), ncp,
                                  ndir, stepSize, tolerance);
  BOOST_CHECK_NE(esState.jacToGlobal, BoundToFreeMatrix::Zero());
  BOOST_CHECK(esState.covTransport);
  BOOST_CHECK_EQUAL(esState.cov, cov);*/
}
