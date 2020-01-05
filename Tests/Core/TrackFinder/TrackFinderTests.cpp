// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE TrackFinder Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/TrackFinder/TrackFinder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

struct ExtendedMinimalSourceLink {
  size_t sourceID;

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

// A few initialisations and definitionas
using SourceLink = ExtendedMinimalSourceLink;
using Jacobian = BoundParameters::CovMatrix_t;
using Covariance = BoundSymMatrix;

using TrackState = TrackState<SourceLink, BoundParameters>;
using Resolution = std::pair<ParID_t, double>;
using ElementResolution = std::vector<Resolution>;
using VolumeResolution = std::map<GeometryID::Value, ElementResolution>;
using DetectorResolution = std::map<GeometryID::Value, VolumeResolution>;

using DebugOutput = detail::DebugOutputActor;

std::normal_distribution<double> gauss(0., 1.);
std::default_random_engine generator(42);

ActsSymMatrixD<1> cov1D;
ActsSymMatrixD<2> cov2D;

bool debugMode = false;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();
CalibrationContext calContext = CalibrationContext();

template <ParID_t... params>
using MeasurementType = Measurement<SourceLink, params...>;

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
      auto geoID = surface->geoID();
      auto volumeID = geoID.volume();
      auto layerID = geoID.layer();
      // find volume and layer information for this
      auto vResolution = detectorResolution.find(volumeID);
      if (vResolution != detectorResolution.end()) {
        // find layer resolutions
        auto lResolution = vResolution->second.find(layerID);
        if (lResolution != vResolution->second.end()) {
          // Apply global to local
          Acts::Vector2D lPos;
          surface->globalToLocal(state.geoContext,
                                 stepper.position(state.stepping),
                                 stepper.direction(state.stepping), lPos);
          if (lResolution->second.size() == 1) {
            double sp = lResolution->second[0].second;
            cov1D << sp * sp;
            double dp = sp * gauss(generator);
            if (lResolution->second[0].first == eLOC_0) {
              // push back & move a LOC_0 measurement
              MeasurementType<eLOC_0> m0(surface->getSharedPtr(), {}, cov1D,
                                         lPos[eLOC_0] + dp);
              result.push_back(std::move(m0));
            } else {
              // push back & move a LOC_1 measurement
              MeasurementType<eLOC_1> m1(surface->getSharedPtr(), {}, cov1D,
                                         lPos[eLOC_1] + dp);
              result.push_back(std::move(m1));
            }
          } else if (lResolution->second.size() == 2) {
            // Create the measurment and move it
            double sx = lResolution->second[eLOC_0].second;
            double sy = lResolution->second[eLOC_1].second;
            cov2D << sx * sx, 0., 0., sy * sy;
            double dx = sx * gauss(generator);
            double dy = sy * gauss(generator);
            // push back & move a LOC_0, LOC_1 measurement
            MeasurementType<eLOC_0, eLOC_1> m01(surface->getSharedPtr(), {},
                                                cov2D, lPos[eLOC_0] + dx,
                                                lPos[eLOC_1] + dy);
            result.push_back(std::move(m01));
          }
        }
      }
    }
  }
};

///
/// @brief Unit test for Kalman fitter with measurements along the x-axis
///
BOOST_AUTO_TEST_CASE(track_finder_zero_field) {
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
  using MeasurementActions = ActionList<MeasurementCreator, DebugOutput>;
  using MeasurementAborters = AbortList<detail::EndOfWorldReached>;

  auto pixelResX = Resolution(eLOC_0, 25_um);
  auto pixelResY = Resolution(eLOC_1, 50_um);
  auto stripResX = Resolution(eLOC_0, 100_um);
  auto stripResY = Resolution(eLOC_1, 150_um);

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
      tgContext, mfContext);
  mOptions.debug = debugMode;
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
  Vector3D mMom(1_GeV, 0., 0);
  for (const auto& [trackID, mPos] : startingPos) {
    SingleCurvilinearTrackParameters<NeutralPolicy> mStart(std::nullopt, mPos,
                                                           mMom, 42_ns);
    // Launch and collect - the measurements
    auto mResult = mPropagator.propagate(mStart, mOptions).value();
    if (debugMode) {
      const auto debugString =
          mResult.template get<DebugOutput::result_type>().debugString;
      std::cout << ">>>> Measurement creation: " << std::endl;
      std::cout << debugString;
    }

    // Extract measurements from result of propagation.
    auto measurementsCreated =
        std::move(mResult.template get<MeasurementCreator::result_type>());
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

  // The TrackFinder - we use the eigen stepper for covariance transport
  // Build navigator for the measurement creatoin
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

  using Updater = GainMatrixUpdater<BoundParameters>;
  using Smoother = GainMatrixSmoother<BoundParameters>;
  using TrackFinder = TrackFinder<RecoPropagator, Updater, Smoother>;

  TrackFinder tFinder(rPropagator,
                      getDefaultLogger("TrackFinder", Logging::VERBOSE));

  // Run the KamanFitter for track finding from different starting parameter
  for (const auto& [trackID, pos] : startingPos) {
    // Set initial parameters for the particle track
    Covariance cov;
    cov << pow(10_um, 2), 0., 0., 0., 0., 0., 0., pow(10_um, 2), 0., 0., 0., 0.,
        0., 0., pow(0.0002, 2), 0., 0., 0., 0., 0., 0., pow(0.0002, 2), 0., 0.,
        0., 0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 1.;

    Vector3D rPos =
        pos + Vector3D{0, 10_um * gauss(generator), 10_um * gauss(generator)};
    double rTheta = 0.0002 * gauss(generator);
    double rPhi = 0.0002 * gauss(generator);
    Vector3D rMom(1_GeV * cos(rTheta) * cos(rPhi),
                  1_GeV * cos(rTheta) * sin(rPhi), 1_GeV * sin(rTheta));

    SingleCurvilinearTrackParameters<ChargedPolicy> rStart(cov, rPos, rMom, 1.,
                                                           42.);

    const Surface* rSurface = &rStart.referenceSurface();

    TrackFinderOptions kfOptions(tgContext, mfContext, calContext, rSurface);

    // Found the track(s)
    auto fitRes = tFinder.findTracks(sourcelinks, rStart, kfOptions);
    BOOST_CHECK(fitRes.ok());
    auto foundTrack = *fitRes;
    auto& fittedStates = foundTrack.fittedStates;
    auto& trackTips = foundTrack.trackTips;

    std::cout << "There are " << trackTips.size()
              << " trajectories found for truth track " << trackID << " : "
              << std::endl;
    size_t iTraj = 0;
    for (const auto& tip : trackTips) {
      std::vector<size_t> sourceIds;
      fittedStates.visitBackwards(tip, [&](const auto& trackState) {
        sourceIds.push_back(trackState.uncalibrated().sourceID);
      });

      BOOST_CHECK_EQUAL(sourceIds.size(), 6);

      std::cout << "The source track id for hits on " << iTraj << " trajectory "
                << " are: " << std::endl;
      size_t numFakeHit = 0;
      for (const auto& id : sourceIds) {
        std::cout << id << " : ";
        numFakeHit = numFakeHit + (id != trackID ? 1 : 0);
      }
      std::cout << std::endl;
      std::cout << "There are " << numFakeHit << " fake hits from other tracks."
                << std::endl;
    }
  }
}

}  // namespace Test
}  // namespace Acts
