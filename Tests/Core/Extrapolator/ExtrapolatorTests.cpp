// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Extrapolator Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/MaterialInteractor.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ExtrapolatorTestGeometry.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Global definitions
  // The path limit abort
  using path_limit = detail::PathLimitReached;

  std::vector<std::shared_ptr<const Surface>> stepState;
  auto tGeometry = testGeometry<PlaneSurface>(stepState);

  // get the navigator and provide the TrackingGeometry
  Navigator navigator(tGeometry);
  using BFieldType          = ConstantBField;
  using EigenStepperType    = EigenStepper<BFieldType>;
  using EigenPropagatorType = Propagator<EigenStepperType, Navigator>;

  const double        Bz = 2. * units::_T;
  BFieldType          bField(0, 0, Bz);
  EigenStepperType    estepper(bField);
  EigenPropagatorType epropagator(std::move(estepper), std::move(navigator));

  const int ntests    = 100;
  bool      debugMode = false;

  // A plane selector for the SurfaceCollector
  struct PlaneSelector
  {

    bool
    operator()(const Surface& sf) const
    {
      return (sf.type() == Surface::Plane);
    }
  };

  // This test case checks that no segmentation fault appears
  // - simple extrapolation test
  BOOST_DATA_TEST_CASE(
      test_extrapolation_,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {

    double dcharge = -1 + 2 * charge;
    (void)index;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);

    PropagatorOptions<> options;
    options.maxStepSize = 10. * units::_cm;
    options.pathLimit   = 25 * units::_cm;

    BOOST_TEST(
        (epropagator.propagate(start, options).endParameters != nullptr));
  }

  // This test case checks that no segmentation fault appears
  // - this tests the collection of surfaces
  BOOST_DATA_TEST_CASE(
      test_surface_collection_,
      bdata::random((bdata::seed = 10,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 11,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 12,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 13,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {

    double dcharge = -1 + 2 * charge;
    (void)index;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);

    // A PlaneSelector for the SurfaceCollector
    using PlaneCollector = SurfaceCollector<PlaneSelector>;

    PropagatorOptions<ActionList<PlaneCollector>> options;

    options.maxStepSize = 10. * units::_cm;
    options.pathLimit   = 25 * units::_cm;

    const auto& result           = epropagator.propagate(start, options);
    auto        collector_result = result.get<PlaneCollector::result_type>();

    // step through the surfaces and go step by step
    PropagatorOptions<> optionsEmpty;

    optionsEmpty.maxStepSize = 25. * units::_cm;
    // try propagation from start to each surface
    for (const auto& colsf : collector_result.collected) {
      // get the surface
      const auto& csurface = colsf.surface;
      const auto& cresult
          = epropagator.propagate(start, *csurface, optionsEmpty).endParameters;
      bool worked = (cresult != nullptr);
      BOOST_TEST(worked);
    }
  }

  // This test case checks that no segmentation fault appears
  // - this tests the collection of surfaces
  BOOST_DATA_TEST_CASE(
      test_material_interactor_,
      bdata::random((bdata::seed = 20,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 21,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 22,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 23,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {

    double dcharge = -1 + 2 * charge;
    (void)index;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);

    using DebugOutput = detail::DebugOutputActor;

    PropagatorOptions<ActionList<MaterialInteractor, DebugOutput>> options;
    options.debug       = debugMode;
    options.maxStepSize = 25. * units::_cm;
    options.pathLimit   = 25 * units::_cm;

    const auto& result = epropagator.propagate(start, options);
    if (result.endParameters) {
      // test that you actually lost some energy
      BOOST_TEST(result.endParameters->momentum().norm()
                 < start.momentum().norm());
    }

    if (debugMode) {
      const auto& output = result.get<DebugOutput::result_type>();
      std::cout << ">>> Extrapolation output " << std::endl;
      std::cout << output.debugString << std::endl;
    }
  }

  // This test case checks that no segmentation fault appears
  // - this tests the loop protection
  BOOST_DATA_TEST_CASE(
      loop_protection_test,
      bdata::random((bdata::seed = 20,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.1 * units::_GeV,
                                                        0.5 * units::_GeV)))
          ^ bdata::random((bdata::seed = 21,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 22,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 23,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(100),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    double dcharge = -1 + 2 * charge;
    (void)index;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);

    // Action list and abort list
    using DebugOutput = detail::DebugOutputActor;

    PropagatorOptions<ActionList<MaterialInteractor, DebugOutput>> options;
    options.maxStepSize = 25. * units::_cm;
    options.pathLimit   = 1500. * units::_mm;

    const auto& status = epropagator.propagate(start, options);
    // this test assumes state.options.loopFraction = 0.5
    // maximum momentum allowed
    double pmax = units::SI2Nat<units::MOMENTUM>(
        options.pathLimit * bField.getField(pos).norm() / M_PI);
    if (mom.norm() < pmax) {
      BOOST_CHECK(status.pathLength < options.pathLimit);
    } else {
      BOOST_CHECK_CLOSE(status.pathLength, options.pathLimit, 1e-3);
    }
  }

}  // namespace Test
}  // namespace Acts
