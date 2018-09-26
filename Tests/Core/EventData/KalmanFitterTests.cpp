// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include <random>
#include <vector>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "DetectorBuild.hpp"

namespace Acts {
namespace Test {

  using id = unsigned long int;

  ///
  /// @brief This struct collects surfaces which are hit by the propagator and
  /// which carries at least one measurement
  ///
  struct SurfaceCollection
  {
    // Collection of measurements sorted by their surfaces
    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements;

    /// @brief Constructor
    SurfaceCollection() = default;

    using result_type = std::vector<Surface const*>;

    /// @brief Operater that is callable by an ActionList. The function collects
    /// the surfaces
    ///
    /// @tparam propagator_state_t Type of the propagator state
    /// @param [in] state State of the propagator
    /// @param [out] result Vector of matching surfaces
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state, result_type& result) const
    {
      //~ std::cout << state.options.debugString << std::endl;
      if (state.navigation.currentSurface)
        std::cout << state.stepping.position().x() << " "
                  << state.stepping.position().y() << " "
                  << state.stepping.position().z() << "\t"
                  << state.navigation.currentSurface->geoID().toString()
                  << std::endl;
      if (measurements.find(state.navigation.currentSurface)
          != measurements.end()) {
        std::cout << "true" << std::endl;
        result.push_back(state.navigation.currentSurface);
      }
    }
  };

  ///
  /// @brief Selector structure for the SurfaceCollector
  ///
  struct SelectSurfaceWithHit
  {
    // Collection of measurements sorted by their surfaces
    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements;

    /// @brief Constructor
    SelectSurfaceWithHit() = default;

    /// @brief Operator that tests if a surface has a measurements
    ///
    /// @param [in] sur Surface that is tested for measurements
    /// @return Boolean result of the test
    bool
    operator()(Surface const& sur) const
    {
      return (measurements.find(&sur) != measurements.end());
    }
  };

  struct EndOfWorld
  {
    EndOfWorld() = default;

    template <typename propagator_state_t>
    bool
    operator()(propagator_state_t& state) const
    {
      if (std::abs(state.stepping.position().x()) > 3. * units::_m
          || std::abs(state.stepping.position().y()) > 0.5 * units::_m
          || std::abs(state.stepping.position().z()) > 0.5 * units::_m)
        return true;
      return false;
    }
  };
  ///
  /// @brief Unit test for Kalman fitter with measurements along the x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_initialization)
  {
    // Build detector
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

    // Construct measurements
    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements;

    ActsSymMatrixD<2> cov2D;
    cov2D << 1. * units::_mm, 0., 0., 1. * units::_mm;

    Vector3D       pos(-2. * units::_m, 0., 0.);
    Surface const* sur = detector->lowestTrackingVolume(pos)
                             ->associatedLayer(pos)
                             ->surfaceArray()
                             ->at(pos)[0];
    measurements[sur].push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 0, cov2D, 0., 0.));
    pos = {-1. * units::_m, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 1, cov2D, 0., 0.));

    ActsSymMatrixD<1> cov1D;
    cov1D << 1. * units::_mm;

    pos = {1. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_0>(*sur, 2, cov1D, 0.));

    pos = {1. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_1>(*sur, 3, cov1D, 0.));

    pos = {2. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_0>(*sur, 4, cov1D, 0.));

    pos = {2. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_1>(*sur, 5, cov1D, 0.));

    // Build navigator
    Navigator navi(detector);
    navi.resolvePassive   = false;
    navi.resolveMaterial  = false;
    navi.resolveSensitive = true;

    // Use default stepper
    StraightLineStepper sls;
    // Build navigator
    Propagator<StraightLineStepper, Navigator> prop(sls, navi);

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(-3. * units::_m, 0., 0.),
        //~ Vector3D startParams(-1.0006 * units::_m, 0., 499.8 * units::_mm),
        startMom(1. * units::_GeV, 0., 0);

    SingleCurvilinearTrackParameters<NeutralPolicy> sbtp(
        std::move(covPtr), startParams, startMom);

    // Create action list for surface collection
    ActionList<SurfaceCollection, SurfaceCollector<SelectSurfaceWithHit>> aList;
    aList.get<SurfaceCollection>().measurements = measurements;
    aList.get<SurfaceCollector<SelectSurfaceWithHit>>().selector.measurements
        = measurements;
    // Set options for propagator
    Propagator<StraightLineStepper, Navigator>::
        Options<ActionList<SurfaceCollection,
                           SurfaceCollector<SelectSurfaceWithHit>>>
            propOpts;
    propOpts.actionList = aList;

    // Launch and collect
    const auto&                        result = prop.propagate(sbtp, propOpts);
    const std::vector<Surface const*>& surResult
        = result.get<typename SurfaceCollection::result_type>();
    const SurfaceCollector<SelectSurfaceWithHit>::this_result& surResult2
        = result.get<
            typename SurfaceCollector<SelectSurfaceWithHit>::result_type>();

    // Test if results match the number of measurements
    BOOST_TEST(surResult.size() == 6);
    BOOST_TEST(surResult2.collected.size() == 6);

    // Re-configure propagation with B-field
    ConstantBField               bField(Vector3D(0., 0.5 * units::_T, 0.));
    EigenStepper<ConstantBField> es(bField);
    Propagator<EigenStepper<ConstantBField>, Navigator> propB(es, navi);
    covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtpB(
        std::move(covPtr), startParams, startMom, 1.);
    AbortList<EndOfWorld> abortList;
    Propagator<EigenStepper<ConstantBField>, Navigator>::
        Options<ActionList<SurfaceCollection,
                           SurfaceCollector<SelectSurfaceWithHit>>,
                AbortList<EndOfWorld>>
            propOptsB;
    propOptsB.actionList     = aList;
    propOptsB.stopConditions = abortList;
    propOptsB.maxSteps       = 1e6;

    const auto& resultB = propB.propagate(sbtpB, propOptsB);
    const std::vector<Surface const*>& surResultB
        = resultB.get<typename SurfaceCollection::result_type>();
    const SurfaceCollector<SelectSurfaceWithHit>::this_result& surResultB2
        = resultB.get<
            typename SurfaceCollector<SelectSurfaceWithHit>::result_type>();

    for (size_t i = 0; i < surResultB2.collected.size(); i++)
      std::cout << surResultB2.collected[i].position.x() << " "
                << surResultB2.collected[i].position.y() << " "
                << surResultB2.collected[i].position.z() << std::endl
                << std::endl;
    BOOST_TEST(surResultB.size() == 2);
    BOOST_TEST(surResultB2.collected.size() == 2);
  }

  ///
  /// @brief Unit test for Kalman fitter with measurements with noise along the
  /// x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_with_noise)
  {
    // Build detector
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

    // Construct measurements
    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements;

    std::normal_distribution<double> gauss(0.0, 2. * units::_cm);
    std::default_random_engine       generator(42);

    ActsSymMatrixD<2> cov2D;
    double dX = gauss(generator), dY = gauss(generator);
    cov2D << dX * dX, 0., 0., dY * dY;

    Vector3D       pos(-2. * units::_m, 0., 0.);
    Surface const* sur = detector->lowestTrackingVolume(pos)
                             ->associatedLayer(pos)
                             ->surfaceArray()
                             ->at(pos)[0];
    measurements[sur].push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 0, cov2D, dX, dY));

	dX = gauss(generator), dY = gauss(generator);
	cov2D << dX * dX, 0., 0., dY * dY;
    pos = {-1. * units::_m, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(
        Measurement<id, eLOC_0, eLOC_1>(*sur, 1, cov2D, dX, dY));

    ActsSymMatrixD<1> cov1D;
    dX = gauss(generator);
    cov1D << dX * dX;

    pos = {1. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_0>(*sur, 2, cov1D, dX));

    dX = gauss(generator);
    cov1D << dX * dX;
    pos = {1. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_1>(*sur, 3, cov1D, dX));

    dX = gauss(generator);
    cov1D << dX * dX;
    pos = {2. * units::_m - 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_0>(*sur, 4, cov1D, dX));

    dX = gauss(generator);
    cov1D << dX * dX;
    pos = {2. * units::_m + 1. * units::_mm, 0., 0.};
    sur = detector->lowestTrackingVolume(pos)
              ->associatedLayer(pos)
              ->surfaceArray()
              ->at(pos)[0];
    measurements[sur].push_back(Measurement<id, eLOC_1>(*sur, 5, cov1D, dX));

    // Build navigator
    Navigator navi(detector);
    navi.resolvePassive   = false;
    navi.resolveMaterial  = false;
    navi.resolveSensitive = true;

    // Use default stepper
    StraightLineStepper sls;
    // Build navigator
    Propagator<StraightLineStepper, Navigator> prop(sls, navi);

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(-3. * units::_m, 0., 0.),
        startMom(1. * units::_GeV, 0., 0);

    SingleCurvilinearTrackParameters<NeutralPolicy> sbtp(
        std::move(covPtr), startParams, startMom);

    // Create action list for surface collection
    ActionList<SurfaceCollection, SurfaceCollector<SelectSurfaceWithHit>> aList;
    aList.get<SurfaceCollection>().measurements = measurements;
    aList.get<SurfaceCollector<SelectSurfaceWithHit>>().selector.measurements
        = measurements;
    // Set options for propagator
    Propagator<StraightLineStepper, Navigator>::
        Options<ActionList<SurfaceCollection,
                           SurfaceCollector<SelectSurfaceWithHit>>>
            propOpts;
    propOpts.actionList = aList;

    // Launch and collect
    const auto&                        result = prop.propagate(sbtp, propOpts);
    const std::vector<Surface const*>& surResult
        = result.get<typename SurfaceCollection::result_type>();
    const SurfaceCollector<SelectSurfaceWithHit>::this_result& surResult2
        = result.get<
            typename SurfaceCollector<SelectSurfaceWithHit>::result_type>();

    // Test if results match the number of measurements
    BOOST_TEST(surResult.size() == 6);
    BOOST_TEST(surResult2.collected.size() == 6);
  }
}  // namespace Test
}  // namespace Acts
