// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include <algorithm>
#include <math.h>
#include <random>
#include <vector>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdator.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/Common/CubicTrackingGeometry.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryID.hpp"

namespace Acts {
namespace Test {

  // A few initialisations and definitionas
  using Identifier = GeometryID;
  using Jacobian   = ActsMatrixD<5, 5>;

  template <ParID_t... params>
  using MeasuredTrackState
      = MeasuredTrackState<Identifier, BoundParameters, Jacobian, params...>;
  using VariantTrackState
      = VariantTrackState<Identifier, BoundParameters, Jacobian>;
  using Resolution         = std::pair<ParID_t, double>;
  using ElementResolution  = std::vector<Resolution>;
  using VolumeResolution   = std::map<geo_id_value, ElementResolution>;
  using DetectorResolution = std::map<geo_id_value, VolumeResolution>;

  using DebugOutput = detail::DebugOutputActor;

  std::normal_distribution<double> gauss(0., 1.);
  std::default_random_engine       generator(42);

  ActsSymMatrixD<1> cov1D;
  ActsSymMatrixD<2> cov2D;

  bool debugMode = false;

  /// @brief This struct creates FittableMeasurements on the
  /// detector surfaces, according to the given smearing xxparameters
  ///
  struct MeasurementCreator
  {
    /// @brief Constructor
    MeasurementCreator() = default;

    /// The detector resolution
    DetectorResolution detectorResolution;

    using result_type = std::vector<VariantTrackState>;

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
      // monitor the current surface
      auto surface = state.navigation.currentSurface;
      if (surface and surface->associatedDetectorElement()) {
        auto         geoID    = surface->geoID();
        geo_id_value volumeID = geoID.value(GeometryID::volume_mask);
        geo_id_value layerID  = geoID.value(GeometryID::layer_mask);
        // find volume and layer information for this
        auto vResolution = detectorResolution.find(volumeID);
        if (vResolution != detectorResolution.end()) {
          // find layer resolutions
          auto lResolution = vResolution->second.find(layerID);
          if (lResolution != vResolution->second.end()) {
            // Apply global to local
            Acts::Vector2D lPos;
            surface->globalToLocal(
                state.stepping.position(), state.stepping.direction(), lPos);
            if (lResolution->second.size() == 1) {
              double sp = lResolution->second[0].second;
              cov1D << sp * sp;
              double dp = sp * gauss(generator);
              if (lResolution->second[0].first == eLOC_0) {
                // push back & move a LOC_0 measurement
                Measurement<Identifier, eLOC_0> m0(
                    *surface, geoID, cov1D, lPos[eLOC_0] + dp);
                result.push_back(MeasuredTrackState<eLOC_0>(std::move(m0)));
              } else {
                // push back & move a LOC_1 measurement
                Measurement<Identifier, eLOC_1> m1(
                    *surface, geoID, cov1D, lPos[eLOC_1] + dp);
                result.push_back(MeasuredTrackState<eLOC_1>(std::move(m1)));
              }
            } else if (lResolution->second.size() == 2) {
              // Create the measurment and move it
              double sx = lResolution->second[eLOC_0].second;
              double sy = lResolution->second[eLOC_1].second;
              cov2D << sx * sx, 0., 0., sy * sy;
              double dx = sx * gauss(generator);
              double dy = sy * gauss(generator);
              // push back & move a LOC_0, LOC_1 measurement
              Measurement<Identifier, eLOC_0, eLOC_1> m01(
                  *surface, geoID, cov2D, lPos[eLOC_0] + dx, lPos[eLOC_1] + dy);
              result.push_back(
                  MeasuredTrackState<eLOC_0, eLOC_1>(std::move(m01)));
            }
          }
        }
      }
    }
  };

  /** Fitter needs that
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
  */
  /** replace with Propagator/detail/EndOfWorld
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
  */

  double         dX, dY;
  Vector3D       pos;
  Surface const* sur;

  ///
  /// @brief Simplified material interaction effect by pure gaussian deflection
  ///
  struct MaterialScattering
  {
    /// @brief Constructor
    MaterialScattering() = default;

    /// @brief Main action list call operator for the scattering on material
    ///
    /// @todo deal momentum in a gaussian way properly
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
            0., 0.017);  //< \approx 1 degree
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
  /// use Propagator to create measurements
  ///
  /// @param [in] detector Detector geometry for surface lookup
  /// @param [in] surfaces Vector of Coordinates referring to surfaces that will
  /// receive measurements
  /// @param [in] dimensions Vector that states if the measurement has a x
  /// and/or y coordinate
  /// @param [in] noise Boolean expression if the measurements receive
  /// underlying noise
  /// @return Map containing the surfaces and the corresponding measurements
  // std::map<Surface const*, std::vector<FittableMeasurement<id>>>
  // createMeasurements(std::shared_ptr<const TrackingGeometry> detector,
  //                   std::vector<Vector3D>&            surfaces,
  //                   std::vector<std::pair<bool, bool>>& dimensions,
  //                   bool noise)
  //{
  //  std::map<Surface const*, std::vector<FittableMeasurement<id>>>
  //  measurements;
  //
  //  // Walk over every surface
  //  for (unsigned long int i = 0; i < surfaces.size(); i++) {
  //    dX = noise ? gauss(generator) : 0.;
  //    // Produce a measurement with x and y coordinate
  //    if (dimensions[i].first && dimensions[i].second) {
  //      dY = noise ? gauss(generator) : 0.;
  //      cov2D << dX * dX, 0., 0., dY * dY;
  //      sur = detector->lowestTrackingVolume(surfaces[i])
  //                ->associatedLayer(surfaces[i])
  //                ->surfaceArray()
  //                ->at(surfaces[i])[0];
  //      measurements[sur].push_back(
  //          Measurement<id, eLOC_0, eLOC_1>(*sur, i, cov2D, dX, dY));
  //    } else {
  //      // Produce measurement with x XOR y coordinate
  //      cov1D << dX * dX;
  //      sur = detector->lowestTrackingVolume(surfaces[i])
  //                ->associatedLayer(surfaces[i])
  //                ->surfaceArray()
  //                ->at(surfaces[i])[0];
  //      if (dimensions[i].first) {
  //        measurements[sur].push_back(
  //            Measurement<id, eLOC_0>(*sur, i, cov1D, dX));
  //      } else if (dimensions[i].second) {
  //        measurements[sur].push_back(
  //            Measurement<id, eLOC_1>(*sur, i, cov1D, dX));
  //      }
  //    }
  //  }
  //
  //  return measurements;
  //}

  ///
  /// @brief Unit test for Kalman fitter with measurements along the x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_zero_field)
  {
    // Build detector
    CubicTrackingGeometry cGeometry;
    auto                  detector = cGeometry();

    // Build navigator for the measurement creatoin
    Navigator mNavigator(detector);
    mNavigator.resolvePassive   = false;
    mNavigator.resolveMaterial  = true;
    mNavigator.resolveSensitive = true;

    // Use straingt line stepper to create the measurements
    StraightLineStepper mStepper;

    // Define the measurement propagator
    using MeasurementPropagator = Propagator<StraightLineStepper, Navigator>;

    // Build propagator for the measurement creation
    MeasurementPropagator mPropagator(mStepper, mNavigator);
    Vector3D mPos(-3. * units::_m, 0., 0.), mMom(1. * units::_GeV, 0., 0);
    SingleCurvilinearTrackParameters<NeutralPolicy> mStart(nullptr, mPos, mMom);

    // Create action list for the measurement creation
    using MeasurementActions  = ActionList<MeasurementCreator, DebugOutput>;
    using MeasurementAborters = AbortList<detail::EndOfWorldReached>;

    auto pixelResX = Resolution(eLOC_0, 25. * units::_um);
    auto pixelResY = Resolution(eLOC_1, 50. * units::_um);
    auto stripResX = Resolution(eLOC_0, 100. * units::_um);
    auto stripResY = Resolution(eLOC_1, 150. * units::_um);

    ElementResolution pixelElementRes  = {pixelResX, pixelResY};
    ElementResolution stripElementResI = {stripResX};
    ElementResolution stripElementResO = {stripResY};

    VolumeResolution pixelVolumeRes;
    pixelVolumeRes[2] = pixelElementRes;
    pixelVolumeRes[4] = pixelElementRes;

    VolumeResolution stripVolumeRes;
    stripVolumeRes[2] = stripElementResI;
    stripVolumeRes[4] = stripElementResO;
    stripVolumeRes[6] = stripElementResI;
    stripVolumeRes[8] = stripElementResO;

    DetectorResolution detRes;
    detRes[2] = pixelVolumeRes;
    detRes[3] = stripVolumeRes;

    // Set options for propagator
    MeasurementPropagator::Options<MeasurementActions, MeasurementAborters>
        mOptions;
    mOptions.debug              = debugMode;
    auto& mCreator              = mOptions.actionList.get<MeasurementCreator>();
    mCreator.detectorResolution = detRes;

    // Launch and collect - the measurements
    auto mResult = mPropagator.propagate(mStart, mOptions);
    if (debugMode) {
      const auto debugString
          = mResult.template get<DebugOutput::result_type>().debugString;
      std::cout << ">>>> Measurement creation: " << std::endl;
      std::cout << debugString;
    }

    auto measurements = mResult.template get<MeasurementCreator::result_type>();
    BOOST_TEST(measurements.size() == 6);

    // The KalmanFitter - we use the eigen stepper for covariance transport
    // Build navigator for the measurement creatoin
    Navigator rNavigator(detector);
    rNavigator.resolvePassive   = false;
    rNavigator.resolveMaterial  = true;
    rNavigator.resolveSensitive = true;

    // Configure propagation with deactivated B-field
    ConstantBField bField(Vector3D(0., 0., 0.));
    using RecoStepper = EigenStepper<ConstantBField>;
    RecoStepper rStepper(bField);
    using RecoPropagator = Propagator<RecoStepper, Navigator>;
    RecoPropagator rPropagator(rStepper, rNavigator);

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 1000. * units::_um, 0., 0., 0., 0., 0., 1000. * units::_um, 0., 0.,
        0., 0., 0., 0.05, 0., 0., 0., 0., 0., 0.05, 0., 0., 0., 0., 0., 0.01;

    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    Vector3D rPos(-3. * units::_m,
                  10. * units::_um * gauss(generator),
                  100. * units::_um * gauss(generator));
    Vector3D rMom(1. * units::_GeV,
                  0.025 * units::_GeV * gauss(generator),
                  0.025 * units::_GeV * gauss(generator));

    SingleCurvilinearTrackParameters<ChargedPolicy> rStart(
        std::move(covPtr), rPos, rMom, 1.);

    const Surface* rSurface = &rStart.referenceSurface();

    using Updator      = GainMatrixUpdator<BoundParameters, Jacobian>;
    using Smoother     = GainMatrixSmoother<BoundParameters, Jacobian>;
    using KalmanFitter = KalmanFitter<RecoPropagator, Updator, Smoother>;

    KalmanFitter kFitter(rPropagator);

    // Fit the track
    auto fittedTrack      = kFitter.fit(measurements, rStart, rSurface);
    auto fittedParameters = fittedTrack.fittedParameters.get();
    /*
    // Make sure it is deterministic
    auto fittedAgainTrack = kFitter.fit(measurements, rStart, rSurface);
    auto fittedAgainParameters = fittedAgainTrack.fittedParameters.get();

    BOOST_TEST(fittedParameters.parameters().isApprox(fittedAgainParameters.parameters()));

    // Change the order of the measurements
    std::vector<VariantTrackState>
      shuffledMeasurements = {  measurements[3],
                                measurements[2],
                                measurements[1],
                                measurements[4],
                                measurements[5],
                                measurements[0] };

    // Make sure it works for shuffled measurements as well
    auto fittedShuffledTrack      = kFitter.fit(shuffledMeasurements, rStart,
   rSurface);
    auto fittedShuffledParameters = fittedShuffledTrack.fittedParameters.get();

    BOOST_TEST(fittedParameters.parameters().isApprox(fittedShuffledParameters.parameters()));

    // Remove one measurement and find a hole
    std::vector<VariantTrackState>
      measurementsWithHole = {  measurements[0],
                                measurements[1],
                                measurements[2],
                                measurements[4],
                                measurements[5] };

   // Make sure it works for shuffled measurements as well
   auto fittedWithHoleTrack      = kFitter.fit(measurementsWithHole, rStart,
   rSurface);
   auto fittedWithHoleParameters = fittedWithHoleTrack.fittedParameters.get();

   // Count one hole
   BOOST_TEST(fittedWithHoleTrack.missedActiveSurfaces.size() == 1);
   // And the parameters should be different
   BOOST_TEST(!fittedParameters.parameters().isApprox(fittedWithHoleParameters.parameters()));
   */
  }

  /*
  ///
  /// @brief Unit test for Kalman fitter with measurements with noise along the
  /// x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_noisy)
  {
    // Build detector
    CubicTrackingGeometry cGeometry;
    auto detector = cGeometry();

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
    const auto& result = prop.propagate(sbtp, propOpts);
    const std::vector<FittableMeasurement<id>>& surResult
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
    const auto& resultB = prop.propagate(sbtp, propOpts);
    const std::vector<FittableMeasurement<id>>& surResultB
        = resultB.get<typename SurfaceCollection::result_type>();
    const SurfaceCollector<SelectSurfaceWithHit>::this_result& surResultB2
        = resultB.get<
            typename SurfaceCollector<SelectSurfaceWithHit>::result_type>();

    BOOST_TEST(surResultB.size() == 4);
    BOOST_TEST(surResultB2.collected.size() == 4);
  }
  */
}  // namespace Test
}  // namespace Acts
