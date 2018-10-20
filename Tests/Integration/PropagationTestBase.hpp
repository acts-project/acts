// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Propagator Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include <cmath>

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

/// test consistency of forward-backward propagation
BOOST_DATA_TEST_CASE(
    forward_backward_propagation_,
    bdata::random((bdata::seed = 0,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 1,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 2,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.1, M_PI - 0.1)))
        ^ bdata::random((bdata::seed = 3,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 4,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0, 1. * units::_m)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    plimit,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;

  // foward backward check atlas stepper
  foward_backward(apropagator,
                  pT,
                  phi,
                  theta,
                  dcharge,
                  plimit,
                  index,
                  1e-3,
                  Acts::units::_eV,
                  debug);
  // foward backward check eigen stepper
  foward_backward(epropagator,
                  pT,
                  phi,
                  theta,
                  dcharge,
                  plimit,
                  index,
                  1e-3,
                  Acts::units::_eV,
                  debug);
}

/// test consistency of propagators when approaching a cylinder
BOOST_DATA_TEST_CASE(
    propagation_to_cylinder_,
    bdata::random((bdata::seed = 1010,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 1111,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 1212,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.1, 0.9 * M_PI)))
        ^ bdata::random((bdata::seed = 1313,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 1414,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.5, 1.)))
        ^ bdata::random((bdata::seed = 1515,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 1616,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 1717,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    rfrac,
    rand1,
    rand2,
    rand3,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;
  // just make sure we can reach it
  double r = rfrac * std::abs(Nat2SI<units::MOMENTUM>(pT) / (1. * Bz));
  r        = (r > 2.5 * Acts::units::_m) ? 2.5 * Acts::units::_m : r;

  // check atlas stepper
  auto a_at_cylinder = to_cylinder(apropagator,
                                   pT,
                                   phi,
                                   theta,
                                   dcharge,
                                   r,
                                   rand1,
                                   rand2,
                                   rand3,
                                   covtpr,
                                   debug);
  // check eigen stepper
  auto e_at_cylinder = to_cylinder(epropagator,
                                   pT,
                                   phi,
                                   theta,
                                   dcharge,
                                   r,
                                   rand1,
                                   rand2,
                                   rand3,
                                   covtpr,
                                   debug);
  CHECK_CLOSE_ABS(e_at_cylinder.first, a_at_cylinder.first, 10. * units::_um);
}

/// test consistency of propagators to a plane
BOOST_DATA_TEST_CASE(
    propagation_to_plane_,
    bdata::random((bdata::seed = 0,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 1,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 2,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0., M_PI)))
        ^ bdata::random((bdata::seed = 3,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 4,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.5, 1.)))
        ^ bdata::random((bdata::seed = 5,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 6,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 7,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    pfrac,
    rand1,
    rand2,
    rand3,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;
  // to a plane with the atlas stepper
  auto a_at_plane
      = to_surface<AtlasPropagatorType, PlaneSurface>(apropagator,
                                                      pT,
                                                      phi,
                                                      theta,
                                                      dcharge,
                                                      pfrac * Acts::units::_m,
                                                      rand1,
                                                      rand2,
                                                      rand3,
                                                      true,
                                                      covtpr);
  // to a plane with the eigen stepper
  auto e_at_plane
      = to_surface<EigenPropagatorType, PlaneSurface>(epropagator,
                                                      pT,
                                                      phi,
                                                      theta,
                                                      dcharge,
                                                      pfrac * Acts::units::_m,
                                                      rand1,
                                                      rand2,
                                                      rand3,
                                                      true,
                                                      covtpr);

  CHECK_CLOSE_ABS(e_at_plane.first, a_at_plane.first, 1 * units::_um);
}

/// test consistency of propagators to a disc
BOOST_DATA_TEST_CASE(
    propagation_to_disc_,
    bdata::random((bdata::seed = 0,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 1,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 2,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.1, M_PI - 0.1)))
        ^ bdata::random((bdata::seed = 3,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 4,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.5, 1.)))
        ^ bdata::random((bdata::seed = 5,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 6,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 7,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    pfrac,
    rand1,
    rand2,
    rand3,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;
  // to a disc with the  atlas stepper
  auto a_at_disc
      = to_surface<AtlasPropagatorType, DiscSurface>(apropagator,
                                                     pT,
                                                     phi,
                                                     theta,
                                                     dcharge,
                                                     pfrac * Acts::units::_m,
                                                     rand1,
                                                     rand2,
                                                     rand3,
                                                     true,
                                                     covtpr);
  // to a disc with the eigen stepper
  auto e_at_disc
      = to_surface<EigenPropagatorType, DiscSurface>(epropagator,
                                                     pT,
                                                     phi,
                                                     theta,
                                                     dcharge,
                                                     pfrac * Acts::units::_m,
                                                     rand1,
                                                     rand2,
                                                     rand3,
                                                     true,
                                                     covtpr);

  CHECK_CLOSE_ABS(e_at_disc.first, a_at_disc.first, 1 * units::_um);
}

/// test consistency of propagators to a line
BOOST_DATA_TEST_CASE(
    propagation_to_line_,
    bdata::random((bdata::seed = 1000,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 1001,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 1002,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.1, M_PI - 0.1)))
        ^ bdata::random((bdata::seed = 1003,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 1004,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.5, 1.)))
        ^ bdata::random((bdata::seed = 1005,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 1006,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 1007,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    pfrac,
    rand1,
    rand2,
    rand3,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;

  // to a line with the atlas stepper
  if (debug) {
    std::cout << "[ >>>> Testing Atlas Propagator <<<< ]" << std::endl;
  }
  auto a_at_line
      = to_surface<AtlasPropagatorType, StrawSurface>(apropagator,
                                                      pT,
                                                      phi,
                                                      theta,
                                                      dcharge,
                                                      pfrac * Acts::units::_m,
                                                      rand1,
                                                      rand2,
                                                      rand3,
                                                      false,
                                                      covtpr,
                                                      debug);
  // to a line with the eigen stepper
  if (debug) {
    std::cout << "[ >>>> Testing Eigen Propagator <<<< ]" << std::endl;
  }
  auto e_at_line
      = to_surface<EigenPropagatorType, StrawSurface>(epropagator,
                                                      pT,
                                                      phi,
                                                      theta,
                                                      dcharge,
                                                      pfrac * Acts::units::_m,
                                                      rand1,
                                                      rand2,
                                                      rand3,
                                                      false,
                                                      covtpr,
                                                      debug);

  CHECK_CLOSE_ABS(e_at_line.first, a_at_line.first, 1 * units::_um);
}

/// test correct covariance transport for curvilinear parameters
/// this test only works within the
/// s_curvilinearProjTolerance (in: Definitions.hpp)
BOOST_DATA_TEST_CASE(
    covariance_transport_curvilinear_curvilinear_,
    bdata::random((bdata::seed = 2000,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 2001,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 2002,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.10, M_PI - 0.10)))
        ^ bdata::random((bdata::seed = 2003,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 2004,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.5, 1.)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    plimit,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;
  // covariance check for eigen stepper
  covariance_curvilinear(
      epropagator, pT, phi, theta, dcharge, plimit * Acts::units::_m, index);
  // covariance check for eigen stepper in dense environment
  // covariance check fo atlas stepper
  covariance_curvilinear(
      apropagator, pT, phi, theta, dcharge, plimit * Acts::units::_m, index);
}

// test correct covariance transport from disc to disc
BOOST_DATA_TEST_CASE(
    covariance_transport_disc_disc_,
    bdata::random((bdata::seed = 3000,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 3001,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 3002,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.1, M_PI - 0.1)))
        ^ bdata::random((bdata::seed = 3003,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 3004,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.5, 1.)))
        ^ bdata::random((bdata::seed = 3005,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 3006,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 3007,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    plimit,
    rand1,
    rand2,
    rand3,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;
  // covariance check for atlas stepper
  covariance_bound<AtlasPropagatorType, DiscSurface, DiscSurface>(
      apropagator,
      pT,
      phi,
      theta,
      dcharge,
      plimit * Acts::units::_m,
      rand1,
      rand2,
      rand3,
      index,
      true,
      true,
      1e-1);

  // covariance check for eigen stepper
  covariance_bound<EigenPropagatorType, DiscSurface, DiscSurface>(
      epropagator,
      pT,
      phi,
      theta,
      dcharge,
      plimit * Acts::units::_m,
      rand1,
      rand2,
      rand3,
      index,
      true,
      true,
      1e-1);
}

// test correct covariance transport from plane to plane
BOOST_DATA_TEST_CASE(
    covariance_transport_plane_plane_,
    bdata::random((bdata::seed = 4000,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 4001,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 4002,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.1, M_PI - 0.1)))
        ^ bdata::random((bdata::seed = 4003,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 4004,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.5, 1.)))
        ^ bdata::random((bdata::seed = 4005,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 4006,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 4007,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    plimit,
    rand1,
    rand2,
    rand3,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;
  // covariance check for atlas stepper
  covariance_bound<AtlasPropagatorType, PlaneSurface, PlaneSurface>(
      apropagator,
      pT,
      phi,
      theta,
      dcharge,
      plimit * Acts::units::_m,
      rand1,
      rand2,
      rand3,
      index);

  // covariance check for eigen stepper
  covariance_bound<EigenPropagatorType, PlaneSurface, PlaneSurface>(
      epropagator,
      pT,
      phi,
      theta,
      dcharge,
      plimit * Acts::units::_m,
      rand1,
      rand2,
      rand3,
      index);
}

// test correct covariance transport from straw to straw
// for straw surfaces the numerical fixture is actually more difficult
// to calculate
BOOST_DATA_TEST_CASE(
    covariance_transport_line_line_,
    bdata::random((bdata::seed = 1000,
                   bdata::distribution
                   = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 1001,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-M_PI, M_PI)))
        ^ bdata::random((bdata::seed = 1002,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.15, M_PI - 0.15)))
        ^ bdata::random((bdata::seed = 1003,
                         bdata::distribution
                         = std::uniform_int_distribution<>(0, 1)))
        ^ bdata::random((bdata::seed = 1004,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.1, 0.2)))
        ^ bdata::random((bdata::seed = 1005,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-0.25, 0.25)))
        ^ bdata::random((bdata::seed = 1006,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-0.25, 0.25)))
        ^ bdata::random((bdata::seed = 1007,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-0.25, 0.25)))
        ^ bdata::xrange(ntests),
    pT,
    phi,
    theta,
    charge,
    plimit,
    rand1,
    rand2,
    rand3,
    index)
{
  if (index < skip) {
    return;
  }

  double dcharge = -1 + 2 * charge;

  // covariance check for atlas stepper
  covariance_bound<AtlasPropagatorType, StrawSurface, StrawSurface>(
      apropagator,
      pT,
      phi,
      theta,
      dcharge,
      plimit * Acts::units::_m,
      rand1,
      rand2,
      rand3,
      index,
      false,
      false,
      1e-1);

  // covariance check for eigen stepper
  covariance_bound<EigenPropagatorType, StrawSurface, StrawSurface>(
      epropagator,
      pT,
      phi,
      theta,
      dcharge,
      plimit * Acts::units::_m,
      rand1,
      rand2,
      rand3,
      index,
      false,
      false,
      1e-1);
}

/// test correct covariance transport for curvilinear parameters in dense
/// environment
/// this test only works within the
/// s_curvilinearProjTolerance (in: Definitions.hpp)
BOOST_DATA_TEST_CASE(
    dense_covariance_transport_curvilinear_curvilinear_,
    bdata::random((bdata::seed = 2000,
                   bdata::distribution
                   = std::uniform_real_distribution<>(3. * units::_GeV,
                                                      10. * units::_GeV)))
        ^ bdata::random((bdata::seed = 2004,
                         bdata::distribution
                         = std::uniform_real_distribution<>(0.5, 1.)))
        ^ bdata::random((bdata::seed = 3005,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 3006,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::random((bdata::seed = 3007,
                         bdata::distribution
                         = std::uniform_real_distribution<>(-1., 1.)))
        ^ bdata::xrange(ntests),
    pT,
    plimit,
    rand1,
    rand2,
    rand3,
    index)
{
  if (index < skip) {
    return;
  }

  // covariance check for eigen stepper in dense environment
  DensePropagatorType dpropagator = setupDensePropagator();
  covariance_curvilinear(
      dpropagator, pT, 0., M_PI / 2., 1, plimit * Acts::units::_m, index);

  covariance_bound<DensePropagatorType, DiscSurface, DiscSurface>(
      dpropagator,
      pT,
      0.,
      M_PI / 2.,
      1,
      plimit * Acts::units::_m,
      rand1,
      rand2,
      rand3,
      index,
      true,
      true,
      1e-1);

  covariance_bound<DensePropagatorType, PlaneSurface, PlaneSurface>(
      dpropagator,
      pT,
      0.,
      M_PI / 2.,
      1,
      plimit * Acts::units::_m,
      rand1,
      rand2,
      rand3,
      index);
}
