// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
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

#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolator/MaterialInteractor.hpp"
#include "ACTS/Extrapolator/Navigator.hpp"
#include "ACTS/Extrapolator/SurfaceCollector.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialCollector.hpp"
#include "ACTS/Propagator/ActionList.hpp"
#include "ACTS/Propagator/EigenStepper.hpp"
#include "ACTS/Propagator/Propagator.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "ExtrapolatorTestGeometry.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Global definitions
  // The path limit abort
  typedef detail::PathLimitReached path_limit;

  typedef ConstantBField                BField_type;
  typedef EigenStepper<BField_type>     EigenStepper_type;
  typedef Propagator<EigenStepper_type> EigenPropagator_type;

  const double         Bz = 0.;  // 2. * units::_T;
  BField_type          bField(0, 0, Bz);
  EigenStepper_type    estepper(bField);
  EigenPropagator_type epropagator(std::move(estepper));

  std::vector<std::unique_ptr<const Surface>> stepState;
  auto tGeometry = testGeometry<PlaneSurface>(stepState);

  const int ntests     = 1;
  bool      debug_mode = true;

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

    // Action list and abort list
    typedef ActionList<Navigator> ActionList_type;
    typedef AbortList<>           AbortConditions_type;

    typename EigenPropagator_type::template Options<ActionList_type,
                                                    AbortConditions_type>
        navigator_options;
    navigator_options.maxStepSize   = 10. * units::_cm;
    navigator_options.maxPathLength = 25 * units::_cm;

    // get the navigator and provide the TrackingGeometry
    auto& navigator            = navigator_options.actionList.get<Navigator>();
    navigator.trackingGeometry = tGeometry;
    navigator.debug            = debug_mode;
    navigator.collectSensitive = true;
    navigator.collectMaterial  = true;
    navigator.collectPassive   = false;

    BOOST_TEST((epropagator.propagate(start, navigator_options).endParameters
                != nullptr));
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
    typedef SurfaceCollector<PlaneSelector> PlaneCollector;

    // Action list and abort list
    typedef ActionList<Navigator, PlaneCollector> ActionList_type;
    typedef AbortList<> AbortConditions_type;

    typename EigenPropagator_type::template Options<ActionList_type,
                                                    AbortConditions_type>
        navigator_options;
    navigator_options.maxStepSize = 10. * units::_cm;

    navigator_options.maxPathLength = 25 * units::_cm;

    // get the navigator and provide the TrackingGeometry
    auto& navigator            = navigator_options.actionList.get<Navigator>();
    navigator.trackingGeometry = tGeometry;
    navigator.debug            = debug_mode;
    navigator.collectSensitive = true;
    navigator.collectMaterial  = true;
    navigator.collectPassive   = false;

    const auto& result = epropagator.propagate(start, navigator_options);

    auto navigator_result = result.get<Navigator::result_type>();
    auto collector_result = result.get<PlaneCollector::result_type>();

    // step through the surfaces and go step by step
    typedef ActionList<> ActionList_empty;
    typename EigenPropagator_type::template Options<ActionList_empty,
                                                    AbortConditions_type>
        options;
    options.maxStepSize = 25. * units::_cm;
    // try propagation from start to each surface
    for (const auto& colsf : collector_result.collected) {
      // get the surface
      const auto& csurface = colsf.surface;
      const auto& cresult
          = epropagator.propagate(start, *csurface, options).endParameters;
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

    // Action list and abort list
    typedef ActionList<Navigator, MaterialInteractor> ActionList_type;
    typedef AbortList<> AbortConditions_type;

    typename EigenPropagator_type::template Options<ActionList_type,
                                                    AbortConditions_type>
        navigator_options;
    navigator_options.maxStepSize   = 25. * units::_cm;
    navigator_options.maxPathLength = 25 * units::_cm;

    // get the navigator and provide the TrackingGeometry
    auto& navigator            = navigator_options.actionList.get<Navigator>();
    navigator.trackingGeometry = tGeometry;
    navigator.debug            = debug_mode;

    const auto& result
        = epropagator.propagate(start, navigator_options).endParameters;
    if (result) {
      // test that you actually lost some energy
      BOOST_TEST(result->momentum().mag() < start.momentum().mag());
    }
  }

}  // namespace Test
}  // namespace Acts
