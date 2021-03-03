// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
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
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <algorithm>
#include <memory>
#include <random>

namespace {

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using KalmanFitter =
    Acts::KalmanFitter<ConstantFieldPropagator, KalmanUpdater, KalmanSmoother>;

/// Find outliers using plain distance for testing purposes.
///
/// In a real setup, the outlier classification can be much more involved, e.g.
/// by computing the weighted distance/ local chi2 value. Here, the purpose is
/// to test that the basic principle works using simplified, synthetic data.
/// Thus, the simplest possible implementation should do.
struct TestOutlierFinder {
  double distanceMax = std::numeric_limits<double>::max();

  /// Classify a measurement as a valid one or an outlier.
  ///
  /// @tparam track_state_t Type of the track state
  /// @param state The track state to classify
  /// @retval False if the measurement is not an outlier
  /// @retval True if the measurement is an outlier
  template <typename track_state_t>
  bool operator()(const track_state_t& state) const {
    // can't determine an outlier w/o a measurement or predicted parameters
    if (not state.hasCalibrated() or not state.hasPredicted()) {
      return false;
    }
    auto residuals = state.calibrated() - state.projector() * state.predicted();
    auto distance = residuals.norm();
    return (distanceMax <= distance);
  }
};

// Construct a straight-line propagator.
StraightPropagator makeStraightPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator navigator(std::move(geo));
  navigator.resolvePassive = false;
  navigator.resolveMaterial = true;
  navigator.resolveSensitive = true;
  Acts::StraightLineStepper stepper;
  return StraightPropagator(std::move(stepper), std::move(navigator));
}

// Construct a propagator using a constant magnetic field along z.
ConstantFieldPropagator makeConstantFieldPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo, double bz) {
  Acts::Navigator navigator(std::move(geo));
  navigator.resolvePassive = false;
  navigator.resolveMaterial = true;
  navigator.resolveSensitive = true;
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
  ConstantFieldStepper stepper(std::move(field));
  return ConstantFieldPropagator(std::move(stepper), std::move(navigator));
}

// Construct initial track parameters.
Acts::CurvilinearTrackParameters makeParameters() {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(-3_m, 0., 0., 42_ns);
  return Acts::CurvilinearTrackParameters(mPos4, 0_degree, 90_degree, 1_GeV,
                                          1_e, cov);
}

const GeometryContext geoCtx;
const MagneticFieldContext magCtx;
const CalibrationContext calCtx;

// detector geometry
CubicTrackingGeometry geometryStore(geoCtx);
const auto geometry = geometryStore();
// expected number of measurements for the given detector
constexpr size_t nMeasurements = 6u;

// detector resolutions
const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                        {25_um, 50_um}};
const MeasurementResolution resStrip0 = {MeasurementType::eLoc0, {100_um}};
const MeasurementResolution resStrip1 = {MeasurementType::eLoc1, {150_um}};
const MeasurementResolutionMap resolutions = {
    {GeometryIdentifier().setVolume(2), resPixel},
    {GeometryIdentifier().setVolume(3).setLayer(2), resStrip0},
    {GeometryIdentifier().setVolume(3).setLayer(4), resStrip1},
    {GeometryIdentifier().setVolume(3).setLayer(6), resStrip0},
    {GeometryIdentifier().setVolume(3).setLayer(8), resStrip1},
};

// simulation propagator
const auto simPropagator = makeStraightPropagator(geometry);

// reconstruction propagator and fitter
const auto kfLogger = getDefaultLogger("KalmanFilter", Logging::INFO);
const auto kfZeroPropagator = makeConstantFieldPropagator(geometry, 0_T);
const auto kfZero = KalmanFitter(kfZeroPropagator);

std::default_random_engine rng(42);

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFittingKalmanFitter)

BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(), VoidOutlierFinder(),
      LoggerWrapper{*kfLogger}, PropagatorPlainOptions());
  // this is the default option. set anyways for consistency
  kfOptions.referenceSurface = nullptr;

  auto res = kfZero.fit(sourceLinks, start, kfOptions);
  BOOST_REQUIRE(res.ok());

  const auto& val = res.value();
  BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
  BOOST_CHECK(not val.fittedParameters);
  BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
  // check the output status flags
  BOOST_CHECK(val.smoothed);
  BOOST_CHECK(not val.reversed);
  BOOST_CHECK(not val.reset);
  BOOST_CHECK(val.finished);
  BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceForward) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  // initial fitter options configured for backward filtereing mode
  // backward filtering requires a reference surface
  KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(), VoidOutlierFinder(),
      LoggerWrapper{*kfLogger}, PropagatorPlainOptions());
  kfOptions.referenceSurface = &start.referenceSurface();
  // this is the default option. set anyways for consistency
  kfOptions.propagatorPlainOptions.direction = forward;

  // regular smoothing
  {
    kfOptions.reversedFiltering = false;
    auto res = kfZero.fit(sourceLinks, start, kfOptions);
    BOOST_CHECK(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
    BOOST_CHECK(val.fittedParameters);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    // check the output status flags
    BOOST_CHECK(val.smoothed);
    BOOST_CHECK(not val.reversed);
    BOOST_CHECK(not val.reset);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
  }
  // reverse filtering instead of smoothing
  {
    kfOptions.reversedFiltering = true;
    auto res = kfZero.fit(sourceLinks, start, kfOptions);
    BOOST_REQUIRE(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
    BOOST_CHECK(val.fittedParameters);
    // check the output status flags
    BOOST_CHECK(not val.smoothed);
    BOOST_CHECK(val.reversed);
    BOOST_CHECK(not val.reset);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
    // count the number of `smoothed` states
    size_t nSmoothed = 0;
    val.fittedStates.visitBackwards(
        val.trackTip,
        [&nSmoothed](const auto& state) { nSmoothed += state.hasSmoothed(); });
    BOOST_CHECK_EQUAL(nSmoothed, sourceLinks.size());
  }
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceBackward) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  // create a track near the tracker exit for outward->inward filtering
  Vector4 posOuter = start.fourPosition(geoCtx);
  posOuter[ePos0] = 3_m;
  CurvilinearTrackParameters startOuter(posOuter, start.unitDirection(),
                                        start.absoluteMomentum(),
                                        start.charge(), start.covariance());

  KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(), VoidOutlierFinder(),
      LoggerWrapper{*kfLogger}, PropagatorPlainOptions());
  kfOptions.referenceSurface = &startOuter.referenceSurface();
  kfOptions.propagatorPlainOptions.direction = backward;

  // regular smoothing
  {
    kfOptions.reversedFiltering = false;
    auto res = kfZero.fit(sourceLinks, startOuter, kfOptions);
    BOOST_CHECK(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
    BOOST_CHECK(val.fittedParameters);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    // check the output status flags
    BOOST_CHECK(val.smoothed);
    BOOST_CHECK(not val.reversed);
    BOOST_CHECK(not val.reset);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
  }
  // reverse filtering instead of smoothing
  {
    kfOptions.reversedFiltering = true;
    auto res = kfZero.fit(sourceLinks, startOuter, kfOptions);
    BOOST_CHECK(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
    BOOST_CHECK(val.fittedParameters);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    // check the output status flags
    BOOST_CHECK(not val.smoothed);
    BOOST_CHECK(val.reversed);
    BOOST_CHECK(not val.reset);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
    // count the number of `smoothed` states
    size_t nSmoothed = 0;
    val.fittedStates.visitBackwards(
        val.trackTip,
        [&nSmoothed](const auto& state) { nSmoothed += state.hasSmoothed(); });
    BOOST_CHECK_EQUAL(nSmoothed, sourceLinks.size());
  }
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceAtExit) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  // create a boundless target surface near the tracker exit
  Vector3 center(3._m, 0., 0.);
  Vector3 normal(1., 0., 0.);
  auto targetSurface = Surface::makeShared<PlaneSurface>(center, normal);

  KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(), VoidOutlierFinder(),
      LoggerWrapper{*kfLogger}, PropagatorPlainOptions());
  kfOptions.referenceSurface = targetSurface.get();

  auto res = kfZero.fit(sourceLinks, start, kfOptions);
  BOOST_REQUIRE(res.ok());

  const auto& val = res.value();
  BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
  BOOST_CHECK(val.fittedParameters);
  BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
  // check the output status flags
  BOOST_CHECK(val.smoothed);
  BOOST_CHECK(not val.reversed);
  BOOST_CHECK(not val.reset);
  BOOST_CHECK(val.finished);
  BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
}

BOOST_AUTO_TEST_CASE(ZeroFieldShuffled) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(), VoidOutlierFinder(),
      LoggerWrapper{*kfLogger}, PropagatorPlainOptions());
  kfOptions.referenceSurface = &start.referenceSurface();

  BoundVector parameters = BoundVector::Zero();

  // fit w/ all hits in order
  {
    auto res = kfZero.fit(sourceLinks, start, kfOptions);
    BOOST_REQUIRE(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
    BOOST_REQUIRE(val.fittedParameters);
    parameters = val.fittedParameters->parameters();
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    // check the output status flags
    BOOST_CHECK(val.smoothed);
    BOOST_CHECK(not val.reversed);
    BOOST_CHECK(not val.reset);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
  }
  // fit w/ all hits in random order
  {
    auto shuffledSourceLinks = sourceLinks;
    std::shuffle(shuffledSourceLinks.begin(), shuffledSourceLinks.end(), rng);
    auto res = kfZero.fit(shuffledSourceLinks, start, kfOptions);
    BOOST_REQUIRE(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
    BOOST_REQUIRE(val.fittedParameters);
    // check consistency w/ un-shuffled measurements
    CHECK_CLOSE_ABS(val.fittedParameters->parameters(), parameters, 1e-5);
    BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
    // check the output status flags
    BOOST_CHECK(val.smoothed);
    BOOST_CHECK(not val.reversed);
    BOOST_CHECK(not val.reset);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
  }
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithHole) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  // fitter options w/o target surface
  KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(), VoidOutlierFinder(),
      LoggerWrapper{*kfLogger}, PropagatorPlainOptions());

  // always keep the first and last measurement. leaving those in seems to not
  // count the respective surfaces as holes.
  for (size_t i = 1u; (i + 1u) < sourceLinks.size(); ++i) {
    // remove the i-th measurement
    auto withHole = sourceLinks;
    withHole.erase(std::next(withHole.begin(), i));
    BOOST_REQUIRE_EQUAL(withHole.size() + 1u, sourceLinks.size());
    BOOST_TEST_INFO("Removed measurement " << i);

    auto res = kfZero.fit(withHole, start, kfOptions);
    BOOST_REQUIRE(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
    BOOST_CHECK(not val.fittedParameters);
    BOOST_CHECK_EQUAL(val.measurementStates, withHole.size());
    // check the output status flags
    BOOST_CHECK(val.smoothed);
    BOOST_CHECK(not val.reversed);
    BOOST_CHECK(not val.reset);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 1u);
  }
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithOutliers) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  const auto& outlierSourceLinks = measurements.outlierSourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);
  BOOST_REQUIRE_EQUAL(outlierSourceLinks.size(), nMeasurements);

  // fitter options w/o target surface. outlier distance is set to be below the
  // default outlier distance in the `MeasurementsCreator`
  KalmanFitterOptions<TestSourceLinkCalibrator, TestOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(),
      TestOutlierFinder{5_mm}, LoggerWrapper{*kfLogger},
      PropagatorPlainOptions());

  for (size_t i = 0; i < sourceLinks.size(); ++i) {
    // replace the i-th measurement with an outlier
    auto withOutlier = sourceLinks;
    withOutlier[i] = outlierSourceLinks[i];
    BOOST_REQUIRE_EQUAL(withOutlier.size(), sourceLinks.size());
    BOOST_TEST_INFO("Replaced measurement " << i << " with outlier");

    auto res = kfZero.fit(withOutlier, start, kfOptions);
    BOOST_REQUIRE(res.ok());

    const auto& val = res.value();
    BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
    // count the number of outliers
    size_t nOutliers = 0;
    val.fittedStates.visitBackwards(
        val.trackTip, [&nOutliers](const auto& state) {
          nOutliers += state.typeFlags().test(TrackStateFlag::OutlierFlag);
        });
    BOOST_CHECK_EQUAL(nOutliers, 1u);
    BOOST_CHECK(not val.fittedParameters);
    BOOST_CHECK_EQUAL(val.measurementStates, withOutlier.size() - 1u);
    // check the output status flags
    BOOST_CHECK(val.smoothed);
    BOOST_CHECK(not val.reversed);
    BOOST_CHECK(not val.reset);
    BOOST_CHECK(val.finished);
    BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
  }
}

// TODO this is not really Kalman fitter specific. is probably better tested
// with a synthetic trajectory.
BOOST_AUTO_TEST_CASE(GlobalCovariance) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  // fitter options w/o target surface
  KalmanFitterOptions<TestSourceLinkCalibrator, VoidOutlierFinder> kfOptions(
      geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(), VoidOutlierFinder(),
      LoggerWrapper{*kfLogger}, PropagatorPlainOptions());

  auto res = kfZero.fit(sourceLinks, start, kfOptions);
  BOOST_REQUIRE(res.ok());

  // Calculate global track parameters covariance matrix
  const auto& val = res.value();
  auto [trackParamsCov, stateRowIndices] =
      detail::globalTrackParametersCovariance(val.fittedStates, val.trackTip);
  BOOST_CHECK_EQUAL(trackParamsCov.rows(), sourceLinks.size() * eBoundSize);
  BOOST_CHECK_EQUAL(stateRowIndices.size(), sourceLinks.size());
  // Each smoothed track state will have eBoundSize rows/cols in the global
  // covariance. stateRowIndices is a map of the starting row/index with the
  // state tip as the key. Thus, the last track state (i.e. the state
  // corresponding val.trackTip) has a starting row/index = eBoundSize *
  // (nMeasurements - 1), i.e. 6*(6-1) = 30.
  BOOST_CHECK_EQUAL(stateRowIndices.at(val.trackTip),
                    eBoundSize * (nMeasurements - 1));
}

BOOST_AUTO_TEST_SUITE_END()
