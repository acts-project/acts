// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include <math.h>
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
      if (measurements.find(state.navigation.currentSurface)
          != measurements.end()) {
        result.push_back(state.navigation.currentSurface);
        std::cout << "Selected surface "
                  << state.navigation.currentSurface->geoID().toString()
                  << " at position (" << state.stepping.position().x() << ", "
                  << state.stepping.position().y() << ", "
                  << state.stepping.position().z() << ")" << std::endl;
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

  ///
  /// @brief Aborter for the case that a particle leaves the detector
  ///
  struct EndOfWorld
  {
    /// @brief Constructor
    EndOfWorld() = default;

    /// @brief Main call operator for the abort operation
    ///
    /// @tparam propagator_state_t State of the propagator
    /// @param [in] state State of the propagation
    /// @return Boolean statement if the particle is still in the detector
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

  std::normal_distribution<double> gauss(0., 2. * units::_cm);
  std::default_random_engine       generator(42);
  ActsSymMatrixD<1>                cov1D;
  ActsSymMatrixD<2>                cov2D;
  double                           dX, dY;
  Vector3D                         pos;
  Surface const*                   sur;

  ///
  /// @brief Simplified material interaction effect by pure gaussian deflection
  ///
  struct MaterialScattering
  {
    /// @brief Constructor
    MaterialScattering() = default;

    /// @brief Main action list call operator for the scattering on material
    ///
    /// @tparam propagator_state_t State of the propagator
    /// @param [in] state State of the propagation
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state) const
    {
      // Check if there is a surface with material and a covariance is existing
      if (state.navigation.currentSurface
          && state.navigation.currentSurface->associatedMaterial()
          && state.stepping.cov != ActsSymMatrixD<5>::Zero()) {
        // Sample angles
        std::normal_distribution<double> scatterAngle(
            0., 0.17);  //< \approx 10 degree
        double dPhi = scatterAngle(generator), dTheta = scatterAngle(generator);

        // Update the covariance
        state.stepping.cov(ePHI, ePHI) += dPhi * dPhi;
        state.stepping.cov(eTHETA, eTHETA) += dTheta * dTheta;

        // Update the angles
        double theta = std::acos(state.stepping.direction().z());
        double phi   = std::atan2(state.stepping.direction().y(),
                                state.stepping.direction().x());

        state.stepping.update(
            state.stepping.position(),
            {std::sin(theta + dTheta) * std::cos(phi + dPhi),
             std::sin(theta + dTheta) * std::sin(phi + dPhi),
             std::cos(theta + dTheta)},
            std::max(state.stepping.p
                         - std::abs(gauss(generator)) * units::_MeV,
                     0.));
      }
    }
  };

  /// @brief Function to calculate measurements with x and/or y coordinates
  ///
  /// @param [in] detector Detector geometry for surface lookup
  /// @param [in] surfaces Vector of Coordinates referring to surfaces that will
  /// receive measurements
  /// @param [in] dimensions Vector that states if the measurement has a x
  /// and/or y coordinate
  /// @param [in] noise Boolean expression if the measurements receive
  /// underlying noise
  /// @return Map containing the surfaces and the corresponding measurements
  std::map<Surface const*, std::vector<FittableMeasurement<id>>>
  createMeasurements(std::shared_ptr<TrackingGeometry> detector,
                     std::vector<Vector3D>&            surfaces,
                     std::vector<std::pair<bool, bool>>& dimensions,
                     bool noise)
  {
    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements;

    // Walk over every surface
    for (unsigned long int i = 0; i < surfaces.size(); i++) {
      dX = noise ? gauss(generator) : 0.;
      // Produce a measurement with x and y coordinate
      if (dimensions[i].first && dimensions[i].second) {
        dY = noise ? gauss(generator) : 0.;
        cov2D << dX * dX, 0., 0., dY * dY;
        sur = detector->lowestTrackingVolume(surfaces[i])
                  ->associatedLayer(surfaces[i])
                  ->surfaceArray()
                  ->at(surfaces[i])[0];
        measurements[sur].push_back(
            Measurement<id, eLOC_0, eLOC_1>(*sur, i, cov2D, dX, dY));
      } else {
        // Produce measurement with x XOR y coordinate
        cov1D << dX * dX;
        sur = detector->lowestTrackingVolume(surfaces[i])
                  ->associatedLayer(surfaces[i])
                  ->surfaceArray()
                  ->at(surfaces[i])[0];
        if (dimensions[i].first) {
          measurements[sur].push_back(
              Measurement<id, eLOC_0>(*sur, i, cov1D, dX));
        } else if (dimensions[i].second) {
          measurements[sur].push_back(
              Measurement<id, eLOC_1>(*sur, i, cov1D, dX));
        }
      }
    }

    return measurements;
  }

  ///
  /// @brief Unit test for Kalman fitter with measurements along the x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_initialization)
  {
    // Build detector
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

    // Construct measurements
    std::vector<Vector3D> surfaces;
    surfaces.push_back({-2. * units::_m, 0., 0.});
    surfaces.push_back({-1. * units::_m, 0., 0.});
    surfaces.push_back({1. * units::_m - 1. * units::_mm, 0., 0.});
    surfaces.push_back({1. * units::_m + 1. * units::_mm, 0., 0.});
    surfaces.push_back({2. * units::_m - 1. * units::_mm, 0., 0.});
    surfaces.push_back({2. * units::_m + 1. * units::_mm, 0., 0.});

    std::vector<std::pair<bool, bool>> dimensions;
    dimensions.push_back(std::make_pair(true, true));
    dimensions.push_back(std::make_pair(true, true));
    dimensions.push_back(std::make_pair(true, false));
    dimensions.push_back(std::make_pair(false, true));
    dimensions.push_back(std::make_pair(true, false));
    dimensions.push_back(std::make_pair(false, true));

    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements
        = createMeasurements(detector, surfaces, dimensions, false);

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

    BOOST_TEST(surResultB.size() == 2);
    BOOST_TEST(surResultB2.collected.size() == 2);
  }

  ///
  /// @brief Unit test for Kalman fitter with measurements with noise along the
  /// x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_noisy)
  {
    // Build detector
    std::shared_ptr<TrackingGeometry> detector = buildGeometry();

    // Construct measurements
    // Get the position of the sensitive surfaces
    std::vector<Vector3D> surfaces;
    surfaces.push_back({-2. * units::_m, 0., 0.});
    surfaces.push_back({-1. * units::_m, 0., 0.});
    surfaces.push_back({1. * units::_m - 1. * units::_mm, 0., 0.});
    surfaces.push_back({1. * units::_m + 1. * units::_mm, 0., 0.});
    surfaces.push_back({2. * units::_m - 1. * units::_mm, 0., 0.});
    surfaces.push_back({2. * units::_m + 1. * units::_mm, 0., 0.});

    // Define the measured components
    std::vector<std::pair<bool, bool>> dimensions;
    dimensions.push_back(std::make_pair(true, true));
    dimensions.push_back(std::make_pair(true, true));
    dimensions.push_back(std::make_pair(true, false));
    dimensions.push_back(std::make_pair(false, true));
    dimensions.push_back(std::make_pair(true, false));
    dimensions.push_back(std::make_pair(false, true));

    std::map<Surface const*, std::vector<FittableMeasurement<id>>> measurements
        = createMeasurements(detector, surfaces, dimensions, true);

    // Build navigator
    Navigator navi(detector);
    navi.resolvePassive   = false;
    navi.resolveMaterial  = false;
    navi.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(-3. * units::_m, 0., 0.),
        startMom(1. * units::_GeV, 0., 0);

    // Create action list for surface collection
    ActionList<MaterialScattering,
               SurfaceCollection,
               SurfaceCollector<SelectSurfaceWithHit>>
        aList;
    aList.get<SurfaceCollection>().measurements = measurements;
    aList.get<SurfaceCollector<SelectSurfaceWithHit>>().selector.measurements
        = measurements;

    // Configure propagation with deactivated B-field
    ConstantBField               bField(Vector3D(0., 0., 0.));
    EigenStepper<ConstantBField> es(bField);
    Propagator<EigenStepper<ConstantBField>, Navigator> prop(es, navi);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);
    AbortList<EndOfWorld> abortList;
    Propagator<EigenStepper<ConstantBField>, Navigator>::
        Options<ActionList<MaterialScattering,
                           SurfaceCollection,
                           SurfaceCollector<SelectSurfaceWithHit>>,
                AbortList<EndOfWorld>>
            propOpts;
    propOpts.actionList     = aList;
    propOpts.stopConditions = abortList;
    propOpts.maxSteps       = 1e6;

    // Launch and collect results
    const auto&                        result = prop.propagate(sbtp, propOpts);
    const std::vector<Surface const*>& surResult
        = result.get<typename SurfaceCollection::result_type>();
    const SurfaceCollector<SelectSurfaceWithHit>::this_result& surResult2
        = result.get<
            typename SurfaceCollector<SelectSurfaceWithHit>::result_type>();

    BOOST_TEST(surResult.size() == 3);
    BOOST_TEST(surResult2.collected.size() == 3);

    // Rebuild for further stability testing without measurements on every
    // surface
    // Remove a single measurement
    std::map<Surface const*, std::vector<FittableMeasurement<id>>>::iterator it;
    it = measurements.find(detector->lowestTrackingVolume(surfaces[1])
                               ->associatedLayer(surfaces[1])
                               ->surfaceArray()
                               ->at(surfaces[1])[0]);
    measurements.erase(it);

    it = measurements.find(detector->lowestTrackingVolume(surfaces[5])
                               ->associatedLayer(surfaces[5])
                               ->surfaceArray()
                               ->at(surfaces[5])[0]);
    measurements.erase(it);

    // Update measurements in the action list
    aList.get<SurfaceCollection>().measurements = measurements;
    aList.get<SurfaceCollector<SelectSurfaceWithHit>>().selector.measurements
        = measurements;
    propOpts.actionList = aList;

    // Configure propagation with deactivated B-field
    prop   = Propagator<EigenStepper<ConstantBField>, Navigator>(es, navi);
    covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    sbtp   = SingleCurvilinearTrackParameters<ChargedPolicy>(
        std::move(covPtr), startParams, startMom, 1.);

    // Launch and collect results
    const auto&                        resultB = prop.propagate(sbtp, propOpts);
    const std::vector<Surface const*>& surResultB
        = resultB.get<typename SurfaceCollection::result_type>();
    const SurfaceCollector<SelectSurfaceWithHit>::this_result& surResultB2
        = resultB.get<
            typename SurfaceCollector<SelectSurfaceWithHit>::result_type>();

    BOOST_TEST(surResultB.size() == 4);
    BOOST_TEST(surResultB2.collected.size() == 4);
  }

}  // namespace Test
}  // namespace Acts
