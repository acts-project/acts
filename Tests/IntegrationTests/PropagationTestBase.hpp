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

// The constant field test
/// test forward propagation in constant magnetic field
BOOST_DATA_TEST_CASE(
    constant_bfieldforward_propagation_,
    bdata::random((bdata::seed = 0,
                   bdata::distribution =
                       std::uniform_real_distribution<>(0.4_GeV, 10_GeV))) ^
        bdata::random((bdata::seed = 1,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 2,
                       bdata::distribution =
                           std::uniform_real_distribution<>(0.1, M_PI - 0.1))) ^
        bdata::random(
            (bdata::seed = 3,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::random(
            (bdata::seed = 4,
             bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, time, index) {
  if (index < skip) {
    return;
  }

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = -1 + 2 * charge;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  CurvilinearParameters startC(std::nullopt, pos, mom, q, time);

  Vector3D dir = mom.normalized();
  FreeVector pars;
  pars << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  FreeParameters startF(std::nullopt, pars);

  // constant field propagation atlas stepper
  auto aposition = constant_field_propagation<CurvilinearParameters>(
      apropagator, startC, pT, phi, theta, Bz);
  // constant field propagation eigen stepper
  auto epositionCC = constant_field_propagation<CurvilinearParameters>(
      epropagator, startC, pT, phi, theta, Bz);
  auto epositionFC = constant_field_propagation<FreeParameters>(
      epropagator, startC, pT, phi, theta, Bz);
  auto epositionCF = constant_field_propagation<CurvilinearParameters>(
      epropagator, startF, pT, phi, theta, Bz);
  auto epositionFF = constant_field_propagation<FreeParameters>(
      epropagator, startF, pT, phi, theta, Bz);
  // check consistency
  CHECK_CLOSE_REL(epositionCC, aposition, 1e-6);
  CHECK_CLOSE_REL(epositionCC, epositionFC, 1e-6);
  CHECK_CLOSE_REL(epositionCC, epositionCF, 1e-6);
  CHECK_CLOSE_REL(epositionCC, epositionFF, 1e-6);
}

/// test consistency of forward-backward propagation
BOOST_DATA_TEST_CASE(forward_backward_propagation_,
                     ds::trackParameters* ds::propagationLimit, pT, phi, theta,
                     charge, plimit) {
  ++itest;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  CurvilinearParameters startC(std::nullopt, pos, mom, q, time);

  Vector3D dir = mom.normalized();
  FreeVector pars;
  pars << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  FreeParameters startF(std::nullopt, pars);

  // foward backward check atlas stepper
  foward_backward<CurvilinearParameters, CurvilinearParameters>(
      apropagator, plimit, startC, 1_um, 1_eV, debug);
  // foward backward check eigen stepper
  foward_backward<CurvilinearParameters, CurvilinearParameters>(
      epropagator, plimit, startC, 1_um, 1_eV, debug);
  foward_backward<FreeParameters, CurvilinearParameters>(
      epropagator, plimit, startC, 1_um, 1_eV, debug);
  foward_backward<CurvilinearParameters, FreeParameters>(
      epropagator, plimit, startC, 1_um, 1_eV, debug);
  foward_backward<FreeParameters, FreeParameters>(epropagator, plimit, startC,
                                                  1_um, 1_eV, debug);
  foward_backward<CurvilinearParameters, CurvilinearParameters>(
      epropagator, plimit, startF, 1_um, 1_eV, debug);
  foward_backward<FreeParameters, CurvilinearParameters>(
      epropagator, plimit, startF, 1_um, 1_eV, debug);
  foward_backward<CurvilinearParameters, FreeParameters>(
      epropagator, plimit, startF, 1_um, 1_eV, debug);
  foward_backward<FreeParameters, FreeParameters>(epropagator, plimit, startF,
                                                  1_um, 1_eV, debug);
  // foward backward check straight line stepper
  foward_backward<CurvilinearParameters, CurvilinearParameters>(
      spropagator, plimit, startC, 1_um, 1_eV, debug);
  foward_backward<FreeParameters, CurvilinearParameters>(
      spropagator, plimit, startC, 1_um, 1_eV, debug);
  foward_backward<CurvilinearParameters, FreeParameters>(
      spropagator, plimit, startC, 1_um, 1_eV, debug);
  foward_backward<FreeParameters, FreeParameters>(spropagator, plimit, startC,
                                                  1_um, 1_eV, debug);
  foward_backward<CurvilinearParameters, CurvilinearParameters>(
      spropagator, plimit, startF, 1_um, 1_eV, debug);
  foward_backward<FreeParameters, CurvilinearParameters>(
      spropagator, plimit, startF, 1_um, 1_eV, debug);
  foward_backward<CurvilinearParameters, FreeParameters>(
      spropagator, plimit, startF, 1_um, 1_eV, debug);
  foward_backward<FreeParameters, FreeParameters>(spropagator, plimit, startF,
                                                  1_um, 1_eV, debug);
}

/// test consistency of propagators when approaching a cylinder
BOOST_DATA_TEST_CASE(propagation_to_cylinder_,
                     ds::trackParameters* ds::propagationFraction ^
                         ds::threeRandom,
                     pT, phi, theta, charge, pfrac, rand1, rand2, rand3) {
  // just make sure we can reach it
  double r = pfrac * std::abs(pT / Bz);
  r = (r > 2.5_m) ? 2.5_m : r;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);

  std::optional<Covariance> covOpt = std::nullopt;
  if (covtpr) {
    Covariance cov;
    // take some major correlations (off-diagonals)
    // clang-format off
    cov <<
     10_mm, 0, 0.123, 0, 0.5, 0,
     0, 10_mm, 0, 0.162, 0, 0,
     0.123, 0, 0.1, 0, 0, 0,
     0, 0.162, 0, 0.1, 0, 0,
     0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     0, 0, 0, 0, 0, 1_us;
    // clang-format on
    covOpt = cov;
  }
  CurvilinearParameters startC(covOpt, pos, mom, q, time);
  NeutralCurvilinearParameters startN(covOpt, pos, mom, time);

  Vector3D dir = mom.normalized();
  FreeVector parsC, parsN;
  parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  FreeParameters startCF(std::nullopt, parsC);
  NeutralFreeParameters startNF(std::nullopt, parsN);

  // check atlas stepper
  auto a_at_cylinder =
      to_cylinder(apropagator, startC, r, rand1, rand2, rand3, debug);
  // check eigen stepper
  auto e_at_cylinder =
      to_cylinder(epropagator, startC, r, rand1, rand2, rand3, debug);
  auto e_free_at_cylinder =
      to_cylinder(epropagator, startCF, r, rand1, rand2, rand3, debug);
  CHECK_CLOSE_ABS(e_at_cylinder.first, a_at_cylinder.first, 10_um);
  CHECK_CLOSE_ABS(e_at_cylinder.first, e_free_at_cylinder.first, 10_um);

  // check without charge
  auto s_at_cylinder =
      to_cylinder(spropagator, startN, r, rand1, rand2, rand3, debug);
  auto s_free_at_cylinder =
      to_cylinder(spropagator, startNF, r, rand1, rand2, rand3, debug);
  e_at_cylinder =
      to_cylinder(epropagator, startN, r, rand1, rand2, rand3, debug);
  e_free_at_cylinder =
      to_cylinder(epropagator, startNF, r, rand1, rand2, rand3, debug);
  CHECK_CLOSE_ABS(s_at_cylinder.first, s_free_at_cylinder.first, 1_um);
  CHECK_CLOSE_ABS(s_at_cylinder.first, e_at_cylinder.first, 1_um);
  CHECK_CLOSE_ABS(s_at_cylinder.first, e_free_at_cylinder.first, 1_um);
}

/// test consistency of propagators to a plane
BOOST_DATA_TEST_CASE(propagation_to_plane_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);

  std::optional<Covariance> covOpt = std::nullopt;
  if (covtpr) {
    Covariance cov;
    // take some major correlations (off-diagonals)
    // clang-format off
    cov <<
     10_mm, 0, 0.123, 0, 0.5, 0,
     0, 10_mm, 0, 0.162, 0, 0,
     0.123, 0, 0.1, 0, 0, 0,
     0, 0.162, 0, 0.1, 0, 0,
     0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     0, 0, 0, 0, 0, 1_us;
    // clang-format on
    covOpt = cov;
  }
  CurvilinearParameters startC(covOpt, pos, mom, q, time);
  NeutralCurvilinearParameters startN(covOpt, pos, mom, time);

  Vector3D dir = mom.normalized();
  FreeVector parsC, parsN;
  parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  FreeParameters startCF(std::nullopt, parsC);
  NeutralFreeParameters startNF(std::nullopt, parsN);

  // to a plane with the atlas stepper
  auto a_at_plane = to_surface<AtlasPropagatorType, PlaneSurface>(
      apropagator, startC, plimit, rand1, rand2, rand3, true);
  // to a plane with the eigen stepper
  auto e_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      epropagator, startC, plimit, rand1, rand2, rand3, true);
  auto e_free_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      epropagator, startCF, plimit, rand1, rand2, rand3, true);
  CHECK_CLOSE_ABS(e_at_plane.first, a_at_plane.first, 1_um);
  CHECK_CLOSE_ABS(e_at_plane.first, e_free_at_plane.first, 1_um);

  // to a plane with the straight line stepper
  auto s_at_plane = to_surface<StraightPropagatorType, PlaneSurface>(
      spropagator, startN, plimit, rand1, rand2, rand3, true);
  auto s_free_at_plane = to_surface<StraightPropagatorType, PlaneSurface>(
      spropagator, startNF, plimit, rand1, rand2, rand3, true);
  // to a plane with the eigen stepper without charge
  e_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      epropagator, startN, plimit, rand1, rand2, rand3, true);
  e_free_at_plane = to_surface<EigenPropagatorType, PlaneSurface>(
      epropagator, startNF, plimit, rand1, rand2, rand3, true);
  CHECK_CLOSE_ABS(s_at_plane.first, s_free_at_plane.first, 1_um);
  CHECK_CLOSE_ABS(s_at_plane.first, e_at_plane.first, 1_um);
  CHECK_CLOSE_ABS(s_at_plane.first, e_free_at_plane.first, 1_um);
}

/// test consistency of propagators to a disc
BOOST_DATA_TEST_CASE(propagation_to_disc_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);

  std::optional<Covariance> covOpt = std::nullopt;
  if (covtpr) {
    Covariance cov;
    // take some major correlations (off-diagonals)
    // clang-format off
    cov <<
     10_mm, 0, 0.123, 0, 0.5, 0,
     0, 10_mm, 0, 0.162, 0, 0,
     0.123, 0, 0.1, 0, 0, 0,
     0, 0.162, 0, 0.1, 0, 0,
     0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     0, 0, 0, 0, 0, 1_us;
    // clang-format on
    covOpt = cov;
  }
  CurvilinearParameters startC(covOpt, pos, mom, q, time);
  NeutralCurvilinearParameters startN(covOpt, pos, mom, time);

  Vector3D dir = mom.normalized();
  FreeVector parsC, parsN;
  parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  FreeParameters startCF(std::nullopt, parsC);
  NeutralFreeParameters startNF(std::nullopt, parsN);

  // to a disc with the  atlas stepper
  auto a_at_disc = to_surface<AtlasPropagatorType, DiscSurface>(
      apropagator, startC, plimit, rand1, rand2, rand3, true);
  // to a disc with the eigen stepper
  auto e_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      epropagator, startC, plimit, rand1, rand2, rand3, true);
  auto e_free_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      epropagator, startCF, plimit, rand1, rand2, rand3, true);
  CHECK_CLOSE_ABS(a_at_disc.first, e_at_disc.first, 1_um);
  CHECK_CLOSE_ABS(a_at_disc.first, e_free_at_disc.first, 1_um);

  // to a disc with the straight line stepper
  auto s_at_disc = to_surface<StraightPropagatorType, DiscSurface>(
      spropagator, startN, plimit, rand1, rand2, rand3, true);
  auto s_free_at_disc = to_surface<StraightPropagatorType, DiscSurface>(
      spropagator, startNF, plimit, rand1, rand2, rand3, true);
  // to a disc with the eigen stepper without charge
  e_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      epropagator, startN, plimit, rand1, rand2, rand3, true);
  e_free_at_disc = to_surface<EigenPropagatorType, DiscSurface>(
      epropagator, startNF, plimit, rand1, rand2, rand3, true);
  CHECK_CLOSE_ABS(s_at_disc.first, s_free_at_disc.first, 1_um);
  CHECK_CLOSE_ABS(s_at_disc.first, e_at_disc.first, 1_um);
  CHECK_CLOSE_ABS(s_at_disc.first, e_free_at_disc.first, 1_um);
}

/// test consistency of propagators to a line
BOOST_DATA_TEST_CASE(propagation_to_line_,
                     ds::trackParameters* ds::propagationLimit ^
                         ds::threeRandom,
                     pT, phi, theta, charge, plimit, rand1, rand2, rand3) {
  std::cout << "Running (first patch) tests : " << itest << std::endl;
  ++itest;

  // define start parameters
  double x = 1.;
  double y = 0.;
  double z = 0.;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = charge;
  double time = 0.;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);

  std::optional<Covariance> covOpt = std::nullopt;
  if (covtpr) {
    Covariance cov;
    // take some major correlations (off-diagonals)
    // clang-format off
    cov <<
     10_mm, 0, 0.123, 0, 0.5, 0,
     0, 10_mm, 0, 0.162, 0, 0,
     0.123, 0, 0.1, 0, 0, 0,
     0, 0.162, 0, 0.1, 0, 0,
     0.5, 0, 0, 0, 1_e / 10_GeV, 0,
     0, 0, 0, 0, 0, 1_us;
    // clang-format on
    covOpt = cov;
  }
  CurvilinearParameters startC(covOpt, pos, mom, q, time);
  NeutralCurvilinearParameters startN(covOpt, pos, mom, time);

  Vector3D dir = mom.normalized();
  FreeVector parsC, parsN;
  parsC << x, y, z, time, dir.x(), dir.y(), dir.z(), q / mom.norm();
  parsN << x, y, z, time, dir.x(), dir.y(), dir.z(), 1. / mom.norm();
  FreeParameters startCF(std::nullopt, parsC);
  NeutralFreeParameters startNF(std::nullopt, parsN);

  // to a line with the atlas stepper
  if (debug) {
    std::cout << "[ >>>> Testing Atlas Propagator <<<< ]" << std::endl;
  }
  auto a_at_line = to_surface<AtlasPropagatorType, StrawSurface>(
      apropagator, startC, plimit, rand1, rand2, rand3, false, debug);
  // to a line with the eigen stepper
  if (debug) {
    std::cout << "[ >>>> Testing Eigen Propagator <<<< ]" << std::endl;
  }
  auto e_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      epropagator, startC, plimit, rand1, rand2, rand3, false, debug);
  auto e_free_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      epropagator, startCF, plimit, rand1, rand2, rand3, false, debug);
  CHECK_CLOSE_ABS(a_at_line.first, e_at_line.first, 10_um);
  CHECK_CLOSE_ABS(a_at_line.first, e_free_at_line.first, 10_um);

  if (debug) {
    std::cout << "[ >>>> Testing Neutral Propagators <<<< ]" << std::endl;
  }
  // to a straw with the straight line stepper
  auto s_at_line = to_surface<StraightPropagatorType, StrawSurface>(
      spropagator, startN, plimit, rand1, rand2, rand3, false, debug);
  auto s_free_at_line = to_surface<StraightPropagatorType, StrawSurface>(
      spropagator, startNF, plimit, rand1, rand2, rand3, false, debug);
  // to a straw with the eigen stepper without charge
  e_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      epropagator, startN, plimit, rand1, rand2, rand3, false, debug);
  e_free_at_line = to_surface<EigenPropagatorType, StrawSurface>(
      epropagator, startNF, plimit, rand1, rand2, rand3, false, debug);
  CHECK_CLOSE_ABS(s_at_line.first, s_free_at_line.first, 1_um);
  CHECK_CLOSE_ABS(s_at_line.first, e_at_line.first, 1_um);
  CHECK_CLOSE_ABS(s_at_line.first, e_free_at_line.first, 1_um);
}

/// test correct covariance transport for curvilinear parameters
/// this test only works within the
/// s_curvilinearProjTolerance (in: Definitions.hpp)
BOOST_DATA_TEST_CASE(covariance_transport_curvilinear_curvilinear_,
                     ds::trackParameters* ds::propagationLimit, pT, phi, theta,
                     charge, plimit) {
  // covariance check for straight line stepper
  CHECK_CLOSE_COVARIANCE(
      covariance_curvilinear<FreeParameters>(rspropagator, start, plimit),
      covariance_curvilinear<FreeParameters>(spropagator, start, plimit),
      1e-3);
  // covariance check for eigen stepper
  CHECK_CLOSE_COVARIANCE(
      covariance_curvilinear<FreeParameters>(repropagator, start, plimit),
      covariance_curvilinear<FreeParameters>(epropagator, start, plimit),
      1e-3);
      
  //~ ///
  //~ /// Free to Curvilinear Tests
  //~ ///
  
  //~ // covariance check for straight line stepper
  //~ CHECK_CLOSE_COVARIANCE(
      //~ covariance_curvilinear<FreeParameters>(rspropagator, start, plimit),
      //~ covariance_curvilinear<FreeParameters>(spropagator, start, plimit),
      //~ 1e-3);
  //~ // covariance check for eigen stepper
  //~ CHECK_CLOSE_COVARIANCE(
      //~ covariance_curvilinear<FreeParameters>(repropagator, start, plimit),
      //~ covariance_curvilinear<FreeParameters>(epropagator, start, plimit),
      //~ 1e-3);
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
