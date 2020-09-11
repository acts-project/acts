// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/TrackFinding/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::makeVector4;

namespace Acts {
namespace Test {

struct ExtendedMinimalSourceLink {
  size_t sourceID = 0;

  const FittableMeasurement<ExtendedMinimalSourceLink>* meas{nullptr};

  bool operator==(const ExtendedMinimalSourceLink& rhs) const {
    return meas == rhs.meas;
  }

  const Surface& referenceSurface() const {
    return *MeasurementHelpers::getSurface(*meas);
  }

  const FittableMeasurement<ExtendedMinimalSourceLink>& operator*() const {
    return *meas;
  }
};

// helper function to create geometry ids
GeometryIdentifier makeId(int volume = 0, int layer = 0, int sensitive = 0) {
  return GeometryIdentifier().setVolume(volume).setLayer(layer).setSensitive(
      sensitive);
}

// A few initialisations and definitionas
using SourceLink = ExtendedMinimalSourceLink;
using Jacobian = BoundMatrix;
using Covariance = BoundSymMatrix;
using Resolution = std::pair<BoundIndices, double>;
using ElementResolution = std::vector<Resolution>;
using VolumeResolution = std::map<GeometryIdentifier::Value, ElementResolution>;
using DetectorResolution =
    std::map<GeometryIdentifier::Value, VolumeResolution>;

std::normal_distribution<double> gauss(0., 1.);
std::default_random_engine generator(42);

ActsSymMatrixD<1> cov1D;
SymMatrix2D cov2D;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();
CalibrationContext calContext = CalibrationContext();

template <BoundIndices... params>
using MeasurementType = Measurement<SourceLink, BoundIndices, params...>;

/// @brief This struct creates FittableMeasurements on the
/// detector surfaces, according to the given smearing xxparameters
///
struct MeasurementCreator {
  /// @brief Constructor
  MeasurementCreator() = default;

  /// The detector resolution
  DetectorResolution detectorResolution;

  using result_type = std::vector<FittableMeasurement<SourceLink>>;

  /// @brief Operater that is callable by an ActionList. The function collects
  /// the surfaces
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [out] result Vector of matching surfaces
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    // monitor the current surface
    auto surface = state.navigation.currentSurface;
    if (surface and surface->associatedDetectorElement()) {
      auto geoID = surface->geometryId();
      auto volumeID = geoID.volume();
      auto layerID = geoID.layer();
      // find volume and layer information for this
      auto vResolution = detectorResolution.find(volumeID);
      if (vResolution != detectorResolution.end()) {
        // find layer resolutions
        auto lResolution = vResolution->second.find(layerID);
        if (lResolution != vResolution->second.end()) {
          // Apply global to local
          auto lpResult = surface->globalToLocal(
              state.geoContext, stepper.position(state.stepping),
              stepper.direction(state.stepping));
          Acts::Vector2D lPos = lpResult.value();
          if (lResolution->second.size() == 1) {
            double sp = lResolution->second[0].second;
            cov1D << sp * sp;
            double dp = sp * gauss(generator);
            if (lResolution->second[0].first == eBoundLoc0) {
              // push back & move a LOC_0 measurement
              MeasurementType<eBoundLoc0> m0(surface->getSharedPtr(), {}, cov1D,
                                             lPos[eBoundLoc0] + dp);
              result.push_back(std::move(m0));
            } else {
              // push back & move a LOC_1 measurement
              MeasurementType<eBoundLoc1> m1(surface->getSharedPtr(), {}, cov1D,
                                             lPos[eBoundLoc1] + dp);
              result.push_back(std::move(m1));
            }
          } else if (lResolution->second.size() == 2) {
            // Create the measurment and move it
            double sx = lResolution->second[eBoundLoc0].second;
            double sy = lResolution->second[eBoundLoc1].second;
            cov2D << sx * sx, 0., 0., sy * sy;
            double dx = sx * gauss(generator);
            double dy = sy * gauss(generator);
            // push back & move a LOC_0, LOC_1 measurement
            MeasurementType<eBoundLoc0, eBoundLoc1> m01(
                surface->getSharedPtr(), {}, cov2D, lPos[eBoundLoc0] + dx,
                lPos[eBoundLoc1] + dy);
            result.push_back(std::move(m01));
          }
        }
      }
    }
  }
};

///
/// @brief Unit test for CombinatorialKalmanFilter with measurements along the
/// x-axis
///
BOOST_AUTO_TEST_CASE(comb_kalman_filter_zero_field) {
  // Build detector
  CubicTrackingGeometry cGeometry(tgContext);
  auto detector = cGeometry();

  // Build navigator for the measurement creatoin
  Navigator mNavigator(detector);
  mNavigator.resolvePassive = false;
  mNavigator.resolveMaterial = true;
  mNavigator.resolveSensitive = true;

  // Use straingt line stepper to create the measurements
  StraightLineStepper mStepper;

  // Define the measurement propagator
  using MeasurementPropagator = Propagator<StraightLineStepper, Navigator>;

  // Build propagator for the measurement creation
  MeasurementPropagator mPropagator(mStepper, mNavigator);

  // Create action list for the measurement creation
  using MeasurementActions = ActionList<MeasurementCreator>;
  using MeasurementAborters = AbortList<EndOfWorldReached>;

  auto pixelResX = Resolution(eBoundLoc0, 25_um);
  auto pixelResY = Resolution(eBoundLoc1, 50_um);
  auto stripResX = Resolution(eBoundLoc0, 100_um);
  auto stripResY = Resolution(eBoundLoc1, 150_um);

  ElementResolution pixelElementRes = {pixelResX, pixelResY};
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
  PropagatorOptions<MeasurementActions, MeasurementAborters> mOptions(
      tgContext, mfContext, getDummyLogger());
  auto& mCreator = mOptions.actionList.get<MeasurementCreator>();
  mCreator.detectorResolution = detRes;

  // This vector owns the measurements
  std::multimap<size_t, FittableMeasurement<SourceLink>> measurements;

  // Make a vector of source links and further processed for KF inputs
  std::vector<SourceLink> sourcelinks;

  // Set the starting positions for propagation
  double eps = 15_mm;
  std::map<size_t, Vector3D> startingPos;
  startingPos.emplace(0, Vector3D{-3_m, 0., 0.});
  startingPos.emplace(1, Vector3D{-3_m, -1.0 * eps, -1.0 * eps});
  startingPos.emplace(2, Vector3D{-3_m, eps, eps});

  // Run the propagation for a few times such that multiple measurements exist
  // on one surface
  // Set the starting momentum for propagation
  for (const auto& [trackID, mPos] : startingPos) {
    Vector4D pos4 = makeVector4(mPos, 42_ns);
    NeutralCurvilinearTrackParameters mStart(pos4, 0_degree, 90_degree,
                                             1 / 1_GeV);
    // Launch and collect - the measurements
    auto result = mPropagator.propagate(mStart, mOptions);
    BOOST_CHECK(result.ok());

    // Extract measurements from result of propagation.
    auto value = std::move(result.value());
    auto measurementsCreated = value.get<MeasurementCreator::result_type>();
    for (auto& meas : measurementsCreated) {
      measurements.emplace(trackID, std::move(meas));
    }
  }

  // Transform the measurments to sourcelinks
  std::transform(measurements.begin(), measurements.end(),
                 std::back_inserter(sourcelinks), [](const auto& m) {
                   return SourceLink{m.first, &m.second};
                 });

  // There should be 18 source links in total
  BOOST_CHECK_EQUAL(sourcelinks.size(), 18);

  // The CombinatorialKalmanFilter - we use the eigen stepper for covariance
  // transport Build navigator for the measurement creatoin
  Navigator rNavigator(detector);
  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  // Configure propagation with deactivated B-field
  ConstantBField bField(Vector3D(0., 0., 0.));
  using RecoStepper = EigenStepper<ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  using Updater = GainMatrixUpdater;
  using Smoother = GainMatrixSmoother;
  using SourceLinkSelector = CKFSourceLinkSelector;
  using CombinatorialKalmanFilter =
      CombinatorialKalmanFilter<RecoPropagator, Updater, Smoother,
                                SourceLinkSelector>;

  // Implement different chi2/nSourceLinks cutoff at different detector level
  // NB: pixel volumeID = 2, strip volumeID= 3
  SourceLinkSelector::Config sourcelinkSelectorConfig = {
      // global default valies
      {makeId(), {8.0, 10}},
      // pixel layer 2 chi2/nSourceLinks cutoff: 8.0/5
      {makeId(2, 2), {8.0, 5}},
      // pixel layer 4 chi2/nSourceLinks cutoff: 7.0/5
      {makeId(2, 4), {7.0, 5}},
      // pixel volume chi2/nSourceLinks cutoff: 7.0/5
      {makeId(2), {7.0, 5}},
      // strip volume chi2/nSourceLinks cutoff: 8.0/5
      {makeId(3), {8.0, 5}},
  };
  CombinatorialKalmanFilter cKF(rPropagator);

  // Run the CombinaltorialKamanFitter for track finding from different starting
  // parameter
  for (const auto& [trackID, pos] : startingPos) {
    // Set initial parameters for the particle track
    Covariance cov;
    cov << pow(10_um, 2), 0., 0., 0., 0., 0., 0., pow(10_um, 2), 0., 0., 0., 0.,
        0., 0., pow(0.0002, 2), 0., 0., 0., 0., 0., 0., pow(0.0002, 2), 0., 0.,
        0., 0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 1.;
    Vector3D rPos =
        pos + Vector3D{0, 10_um * gauss(generator), 10_um * gauss(generator)};
    double rPhi = 0_degree + 0.0002 * gauss(generator);
    double rTheta = 90_degree + 0.0002 * gauss(generator);
    CurvilinearTrackParameters rStart(makeVector4(rPos, 42_ns), rPhi, rTheta,
                                      1_GeV, 1_e, cov);

    const Surface* rSurface = &rStart.referenceSurface();

    auto logger =
        getDefaultLogger("CombinatorialKalmanFilter", Logging::VERBOSE);
    CombinatorialKalmanFilterOptions<SourceLinkSelector> ckfOptions(
        tgContext, mfContext, calContext, sourcelinkSelectorConfig,
        LoggerWrapper{*logger}, PropagatorPlainOptions(), rSurface);

    // Found the track(s)
    auto combKalmanFilterRes = cKF.findTracks(sourcelinks, rStart, ckfOptions);
    BOOST_CHECK(combKalmanFilterRes.ok());

    auto foundTrack = *combKalmanFilterRes;
    auto& fittedStates = foundTrack.fittedStates;
    auto& trackTips = foundTrack.trackTips;

    for (const auto& tip : trackTips) {
      std::vector<size_t> sourceIds;
      fittedStates.visitBackwards(tip, [&](const auto& trackState) {
        sourceIds.push_back(trackState.uncalibrated().sourceID);
      });

      BOOST_CHECK_EQUAL(sourceIds.size(), 6);

      size_t numFakeHit = 0;
      for (const auto& id : sourceIds) {
        numFakeHit = numFakeHit + (id != trackID ? 1 : 0);
      }

      // Check if there are fake hits from other tracks
      BOOST_CHECK_EQUAL(numFakeHit, 0);
    }
  }
}

}  // namespace Test
}  // namespace Acts
