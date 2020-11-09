// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

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
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/TrackFinding/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <random>
#include <vector>

namespace {

using namespace Acts::Test;
using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using ConstantFieldStepper = Acts::EigenStepper<Acts::ConstantBField>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using CombinatorialKalmanFilter =
    Acts::CombinatorialKalmanFilter<ConstantFieldPropagator, KalmanUpdater,
                                    KalmanSmoother>;

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
  Acts::ConstantBField field(Acts::Vector3D(0.0, 0.0, bz));
  ConstantFieldStepper stepper(std::move(field));
  return ConstantFieldPropagator(std::move(stepper), std::move(navigator));
}

// Construct multiple initial track parameters.
std::vector<Acts::CurvilinearTrackParameters> makeParameters() {
  // create common covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // all tracks close to the transverse plane along the x axis w/ small
  // variations in position, direction.
  Acts::Vector4D mPos0(-3_m, 0.0, 0.0, 1_ns);
  Acts::Vector4D mPos1(-3_m, -15_mm, -15_mm, 2_ns);
  Acts::Vector4D mPos2(-3_m, 15_mm, 15_mm, -1_ns);
  return {
      {mPos0, 0_degree, 90_degree, 1_GeV, 1_e, cov},
      {mPos1, -1_degree, 91_degree, 1_GeV, 1_e, cov},
      {mPos2, 1_degree, 89_degree, 1_GeV, -1_e, cov},
  };
}

const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;
const Acts::CalibrationContext calCtx;

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
    {Acts::GeometryIdentifier().setVolume(2), resPixel},
    {Acts::GeometryIdentifier().setVolume(3).setLayer(2), resStrip0},
    {Acts::GeometryIdentifier().setVolume(3).setLayer(4), resStrip1},
    {Acts::GeometryIdentifier().setVolume(3).setLayer(6), resStrip0},
    {Acts::GeometryIdentifier().setVolume(3).setLayer(8), resStrip1},
};

// source link selection configurtion.
const Acts::CKFSourceLinkSelector::Config selectBestSourceLinkOnly = {
    // global default: no chi2 cut, only one source link per surface
    {Acts::GeometryIdentifier(), {std::numeric_limits<double>::max(), 1u}},
};

// simulation propagator
const auto simPropagator = makeStraightPropagator(geometry);

// reconstruction propagator and fitter
const auto ckfLogger =
    Acts::getDefaultLogger("KalmanFilter", Acts::Logging::INFO);
const auto ckfZeroPropagator = makeConstantFieldPropagator(geometry, 0_T);
const auto ckfZero = CombinatorialKalmanFilter(ckfZeroPropagator);

std::default_random_engine rng(4212352341961);

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFindingCombinatorialKalmanFilter)

BOOST_AUTO_TEST_CASE(ZeroField) {
  auto tracks = makeParameters();

  // simulate all tracks and generate a combined set of resulting sourcelinks
  std::vector<std::shared_ptr<Acts::FittableMeasurement<TestSourceLink>>> store;
  std::vector<TestSourceLink> sourceLinks;
  for (size_t trackId = 0u; trackId < tracks.size(); ++trackId) {
    auto measurements =
        createMeasurements(simPropagator, geoCtx, magCtx, tracks[trackId],
                           resolutions, rng, trackId);
    BOOST_REQUIRE_EQUAL(measurements.sourceLinks.size(), nMeasurements);

    for (auto& fm : measurements.store) {
      store.emplace_back(std::move(fm));
    }
    for (auto& sl : measurements.sourceLinks) {
      sourceLinks.emplace_back(std::move(sl));
    }
  }

  // use one of the parameter surfaces as reference surface for the
  // all reconstructed tracks
  Acts::CombinatorialKalmanFilterOptions<TestSourceLinkCalibrator,
                                         Acts::CKFSourceLinkSelector>
      ckfOptions(geoCtx, magCtx, calCtx, TestSourceLinkCalibrator(),
                 Acts::CKFSourceLinkSelector(selectBestSourceLinkOnly),
                 Acts::LoggerWrapper{*ckfLogger},
                 Acts::PropagatorPlainOptions(),
                 &tracks.front().referenceSurface());

  // run the CKF for each initial track state
  for (size_t trackId = 0u; trackId < tracks.size(); ++trackId) {
    // find the tracks
    auto res = ckfZero.findTracks(sourceLinks, tracks[trackId], ckfOptions);
    BOOST_REQUIRE(res.ok());

    auto val = *res;
    // with the given source link selection cuts, only one trajectory for the
    // given input parameters should be found.
    BOOST_CHECK_EQUAL(val.trackTips.size(), 1u);
    // check purity of first found track
    // find the number of hits not originating from the right track
    size_t numHits = 0u;
    size_t numMissmatchedHits = 0u;
    val.fittedStates.visitBackwards(
        val.trackTips.front(), [&](const auto& trackState) {
          numHits += 1u;
          numMissmatchedHits += (trackId != trackState.uncalibrated().sourceId);
        });
    BOOST_CHECK_EQUAL(numHits, nMeasurements);
    BOOST_CHECK_EQUAL(numMissmatchedHits, 0u);
  }
}

BOOST_AUTO_TEST_SUITE_END()
