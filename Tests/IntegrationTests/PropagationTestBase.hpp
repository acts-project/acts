// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <cmath>

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::UnitLiterals;

unsigned int itest = 0;

// datasets for Boost::Test data-driven test cases
namespace ds {
// track parameters
auto pT = bdata::xrange(0.5_GeV, 10_GeV, 500_MeV);
auto phi = bdata::xrange(-180_degree, 180_degree, 30_degree);
auto theta = bdata::xrange(15_degree, 90_degree, 15_degree);
auto charge = bdata::make({1_e, -1_e});
// combined track parameters as cartesian product of all possible combinations
auto trackParameters = (pT * phi * theta * charge);
auto propagationFraction = bdata::xrange(0.0, 1.0, 0.25);
auto propagationLimit = bdata::xrange(10_cm, 1_m, 10_cm);
// additional random numbers
auto rand1 = bdata::random(
    (bdata::seed = 1515,
     bdata::distribution = std::uniform_real_distribution<>(-1., 1.)));
auto rand2 = bdata::random(
    (bdata::seed = 1616,
     bdata::distribution = std::uniform_real_distribution<>(-1., 1.)));
auto rand3 = bdata::random(
    (bdata::seed = 1717,
     bdata::distribution = std::uniform_real_distribution<>(-1., 1.)));
auto threeRandom = (rand1 ^ rand2 ^ rand2);
}  // namespace ds

/// test consistency of forward-backward propagation
BOOST_DATA_TEST_CASE(forward_backward_propagation_,
                     ds::trackParameters* ds::propagationLimit, pT, phi, theta,
                     charge, plimit) {
  ++itest;

  // foward backward check atlas stepper
  foward_backward(apropagator, pT, phi, theta, charge, plimit, 1_um, 1_eV,
                  debug);
  // foward backward check eigen stepper
  foward_backward(epropagator, pT, phi, theta, charge, plimit, 1_um, 1_eV,
                  debug);
  // foward backward check straight line stepper
  foward_backward(spropagator, pT, phi, theta, charge, plimit, 1_um, 1_eV,
                  debug);
}

/// test consistency of propagators when approaching a cylinder
BOOST_DATA_TEST_CASE(propagation_to_cylinder_,
                     ds::trackParameters* ds::propagationFraction ^
                         ds::threeRandom,
                     pT, phi, theta, charge, pfrac, rand1, rand2, rand3) {
  // just make sure we can reach it
  double r = pfrac * std::abs(pT / Bz);
  r = (r > 2.5_m) ? 2.5_m : r;
  // check atlas stepper
  auto a_at_cylinder = to_cylinder(apropagator, pT, phi, theta, charge, r,
                                   rand1, rand2, rand3, covtpr, debug);
  // check eigen stepper
  auto e_at_cylinder = to_cylinder(epropagator, pT, phi, theta, charge, r,
                                   rand1, rand2, rand3, covtpr, debug);
  CHECK_CLOSE_ABS(e_at_cylinder.first, a_at_cylinder.first, 10_um);

  // check without charge
  auto s_at_cylinder = to_cylinder(spropagator, pT, phi, theta, 0., r, rand1,
                                   rand2, rand3, covtpr, debug);
  e_at_cylinder = to_cylinder(epropagator, pT, phi, theta, 0., r, rand1, rand2,
                              rand3, covtpr, debug);

  CHECK_CLOSE_ABS(s_at_cylinder.first, e_at_cylinder.first, 1_um);
}

/// test consistency of propagators to a plane
BOOST_DATA_TEST_CASE(propagation_to_plane_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  // to a plane with the atlas stepper
  auto a_at_plane = to_surface<AtlasPropagatorType, PlaneSurface>(
      apropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3, true,
      covtpr);
  // to a plane with the eigen stepper
  auto e_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      epropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3, true,
      covtpr);
  CHECK_CLOSE_ABS(e_at_plane.first, a_at_plane.first, 1_um);

  // to a plane with the straight line stepper
  auto s_at_plane = to_surface<StraightPropagatorType, PlaneSurface>(
      spropagator, pT, phi, theta, 0., plimit, rand1, rand2, rand3, true,
      covtpr);
  // to a plane with the eigen stepper without charge
  e_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      epropagator, pT, phi, theta, 0., plimit, rand1, rand2, rand3, true,
      covtpr);
  CHECK_CLOSE_ABS(e_at_plane.first, s_at_plane.first, 1_um);
}

/// test consistency of propagators to a disc
BOOST_DATA_TEST_CASE(propagation_to_disc_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  // to a disc with the  atlas stepper
  auto a_at_disc = to_surface<AtlasPropagatorType, DiscSurface>(
      apropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3, true,
      covtpr);
  // to a disc with the eigen stepper
  auto e_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      epropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3, true,
      covtpr);
  CHECK_CLOSE_ABS(e_at_disc.first, a_at_disc.first, 1_um);

  // to a disc with the straight line stepper
  auto s_at_disc = to_surface<StraightPropagatorType, DiscSurface>(
      spropagator, pT, phi, theta, 0., plimit, rand1, rand2, rand3, true,
      covtpr);
  // to a disc with the eigen stepper without charge
  e_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      epropagator, pT, phi, theta, 0., plimit, rand1, rand2, rand3, true,
      covtpr);

  CHECK_CLOSE_ABS(e_at_disc.first, s_at_disc.first, 1_um);
}

/// test consistency of propagators to a line
BOOST_DATA_TEST_CASE(propagation_to_line_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  // to a line with the atlas stepper
  if (debug) {
    std::cout << "[ >>>> Testing Atlas Propagator <<<< ]" << std::endl;
  }
  auto a_at_line = to_surface<AtlasPropagatorType, StrawSurface>(
      apropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3, false,
      covtpr, debug);
  // to a line with the eigen stepper
  if (debug) {
    std::cout << "[ >>>> Testing Eigen Propagator <<<< ]" << std::endl;
  }
  auto e_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      epropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3, false,
      covtpr, debug);
  CHECK_CLOSE_ABS(e_at_line.first, a_at_line.first, 10_um);

  if (debug) {
    std::cout << "[ >>>> Testing Neutral Propagators <<<< ]" << std::endl;
  }
  // to a straw with the straight line stepper
  auto s_at_line = to_surface<StraightPropagatorType, StrawSurface>(
      spropagator, pT, phi, theta, 0., plimit, rand1, rand2, rand3, false,
      covtpr, debug);
  // to a straw with the eigen stepper without charge
  e_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      epropagator, pT, phi, theta, 0., plimit, rand1, rand2, rand3, false,
      covtpr, debug);

  CHECK_CLOSE_ABS(e_at_line.first, s_at_line.first, 1_um);
}

/// test correct covariance transport for curvilinear parameters
/// this test only works within the
/// s_curvilinearProjTolerance (in: Definitions.hpp)
BOOST_DATA_TEST_CASE(covariance_transport_curvilinear_curvilinear_,
                     ds::trackParameters* ds::propagationLimit, pT, phi, theta,
                     charge, plimit) {
  // covariance check for straight line stepper
  CHECK_CLOSE_COVARIANCE(
      covariance_curvilinear(rspropagator, pT, phi, theta, charge, plimit),
      covariance_curvilinear(spropagator, pT, phi, theta, charge, plimit),
      1e-3);
  // covariance check for eigen stepper
  CHECK_CLOSE_COVARIANCE(
      covariance_curvilinear(repropagator, pT, phi, theta, charge, plimit),
      covariance_curvilinear(epropagator, pT, phi, theta, charge, plimit),
      1e-3);
  // covariance check fo atlas stepper
  CHECK_CLOSE_COVARIANCE(
      covariance_curvilinear(rapropagator, pT, phi, theta, charge, plimit),
      covariance_curvilinear(apropagator, pT, phi, theta, charge, plimit),
      1e-3);
}

// test correct covariance transport from disc to disc
BOOST_DATA_TEST_CASE(covariance_transport_disc_disc_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  // covariance check for straight line stepper
  auto covCalculated =
      covariance_bound<RiddersStraightPropagatorType, DiscSurface, DiscSurface>(
          rspropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          true, true);
  auto covObtained =
      covariance_bound<StraightPropagatorType, DiscSurface, DiscSurface>(
          spropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          true, true);
  if (covCalculated != Covariance::Zero()) {
    CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 5e-1);
  }

  // covariance check for eigen stepper
  covCalculated =
      covariance_bound<RiddersEigenPropagatorType, DiscSurface, DiscSurface>(
          repropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          true, true);
  covObtained = covariance_bound<EigenPropagatorType, DiscSurface, DiscSurface>(
      epropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3, true,
      true);
  if (covCalculated != Covariance::Zero()) {
    CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);
  }

  // covariance check for atlas stepper
  covCalculated =
      covariance_bound<RiddersAtlasPropagatorType, DiscSurface, DiscSurface>(
          rapropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          true, true);
  covObtained = covariance_bound<AtlasPropagatorType, DiscSurface, DiscSurface>(
      apropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3, true,
      true);
  if (covCalculated != Covariance::Zero()) {
    CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);
  }
}

// test correct covariance transport from plane to plane
BOOST_DATA_TEST_CASE(covariance_transport_plane_plane_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  // covariance check for straight line stepper
  auto covCalculated = covariance_bound<RiddersStraightPropagatorType,
                                        PlaneSurface, PlaneSurface>(
      rspropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3);
  auto covObtained =
      covariance_bound<StraightPropagatorType, PlaneSurface, PlaneSurface>(
          spropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3);
  CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-3);

  // covariance check for eigen stepper
  covCalculated =
      covariance_bound<RiddersEigenPropagatorType, PlaneSurface, PlaneSurface>(
          repropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3);
  covObtained =
      covariance_bound<EigenPropagatorType, PlaneSurface, PlaneSurface>(
          epropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3);
  CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-2);

  // covariance check for atlas stepper
  covCalculated =
      covariance_bound<RiddersAtlasPropagatorType, PlaneSurface, PlaneSurface>(
          rapropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3);
  covObtained =
      covariance_bound<AtlasPropagatorType, PlaneSurface, PlaneSurface>(
          apropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3);
  CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-2);
}

// test correct covariance transport from straw to straw
// for straw surfaces the numerical fixture is actually more difficult
// to calculate
BOOST_DATA_TEST_CASE(covariance_transport_line_line_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  // covariance check for straight line stepper
  auto covCalculated =
      covariance_bound<RiddersStraightPropagatorType, StrawSurface,
                       StrawSurface>(rspropagator, pT, phi, theta, charge,
                                     plimit, rand1, rand2, rand3, false, false);
  auto covObtained =
      covariance_bound<StraightPropagatorType, StrawSurface, StrawSurface>(
          spropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          false, false);
  CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);

  // covariance check for eigen stepper
  covCalculated =
      covariance_bound<RiddersEigenPropagatorType, StrawSurface, StrawSurface>(
          repropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          false, false);
  covObtained =
      covariance_bound<EigenPropagatorType, StrawSurface, StrawSurface>(
          epropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          false, false);
  CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);

  // covariance check for atlas stepper
  covCalculated =
      covariance_bound<RiddersAtlasPropagatorType, StrawSurface, StrawSurface>(
          rapropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          false, false);
  covObtained =
      covariance_bound<AtlasPropagatorType, StrawSurface, StrawSurface>(
          apropagator, pT, phi, theta, charge, plimit, rand1, rand2, rand3,
          false, false);
  CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 1e-1);
}

/// test correct covariance transport for curvilinear parameters in dense
/// environment
/// this test only works within the
/// s_curvilinearProjTolerance (in: Definitions.hpp)
BOOST_DATA_TEST_CASE(dense_covariance_transport_curvilinear_curvilinear_,
                     ds::pT* ds::propagationLimit ^ ds::threeRandom, pT, plimit,
                     rand1, rand2, rand3) {
  // covariance check for eigen stepper in dense environment
  DensePropagatorType dpropagator = setupDensePropagator();
  RiddersPropagator rdpropagator(dpropagator);

  CHECK_CLOSE_COVARIANCE(
      covariance_curvilinear(rdpropagator, pT, 0_degree, 45_degree, 1_e,
                             plimit),
      covariance_curvilinear(dpropagator, pT, 0_degree, 45_degree, 1_e, plimit),
      3e-1);

  auto covCalculated =
      covariance_bound<RiddersPropagator<DensePropagatorType>, DiscSurface,
                       DiscSurface>(rdpropagator, pT, 0_degree, 45_degree, 1_e,
                                    plimit, rand1, rand2, rand3, true, true);
  auto covObtained =
      covariance_bound<DensePropagatorType, DiscSurface, DiscSurface>(
          dpropagator, pT, 0_degree, 45_degree, 1_e, plimit, rand1, rand2,
          rand3, true, true);
  if (covCalculated != Covariance::Zero()) {
    CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 8e-1);
  }

  covCalculated = covariance_bound<RiddersPropagator<DensePropagatorType>,
                                   PlaneSurface, PlaneSurface>(
      rdpropagator, pT, 0_degree, 45_degree, 1, plimit, rand1, rand2, rand3);
  covObtained =
      covariance_bound<DensePropagatorType, PlaneSurface, PlaneSurface>(
          dpropagator, pT, 0_degree, 45_degree, 1, plimit, rand1, rand2, rand3);
  CHECK_CLOSE_COVARIANCE(covCalculated, covObtained, 8e-1);
}
