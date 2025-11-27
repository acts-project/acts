// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <functional>
#include <memory>
#include <optional>
#include <random>
#include <utility>

#include "FitterTestsCommon.hpp"

using namespace Acts;
using namespace Acts::detail::Test;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using KalmanFitter =
    Acts::KalmanFitter<ConstantFieldPropagator, VectorMultiTrajectory>;

static const auto pion = Acts::ParticleHypothesis::pion();

KalmanUpdater kfUpdater;

// Construct initial track parameters.
Acts::BoundTrackParameters makeParameters() {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSquareMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(-3_m, 0., 0., 42_ns);
  return Acts::BoundTrackParameters::createCurvilinear(
      mPos4, 0_degree, 90_degree, 1_e / 1_GeV, cov, pion);
}

// Instantiate the tester
const FitterTester tester;

// reconstruction propagator and fitter
auto kfLogger = getDefaultLogger("KalmanFilter", Logging::INFO);
const auto kfZeroPropagator =
    makeConstantFieldPropagator<ConstantFieldStepper>(tester.geometry, 0_T);
const auto kfZero = KalmanFitter(kfZeroPropagator, std::move(kfLogger));

std::default_random_engine rng(42);

auto makeDefaultKalmanFitterOptions() {
  KalmanFitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  extensions.updater.connect<&KalmanUpdater::operator()<VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.surfaceAccessor.connect<
      &Acts::detail::Test::TestSourceLink::SurfaceAccessor::operator()>(
      &tester.surfaceAccessor);

  return KalmanFitterOptions(
      tester.geoCtx, tester.magCtx, tester.calCtx, extensions,
      PropagatorPlainOptions(tester.geoCtx, tester.magCtx));
}

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldNoSurfaceForward(kfZero, kfOptions, start, rng, false,
                                        false, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceForward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldWithSurfaceForward(kfZero, kfOptions, start, rng, false,
                                          false, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceBackward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldWithSurfaceBackward(kfZero, kfOptions, start, rng, false,
                                           false, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceAtExit) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldWithSurfaceAtExit(kfZero, kfOptions, start, rng, false,
                                         false, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldShuffled) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldShuffled(kfZero, kfOptions, start, rng, false, false,
                                true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithHole) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldWithHole(kfZero, kfOptions, start, rng, false, false,
                                true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithOutliers) {
  auto start = makeParameters();

  // fitter options w/o target surface. outlier distance is set to be below the
  // default outlier distance in the `MeasurementsCreator`
  auto kfOptions = makeDefaultKalmanFitterOptions();

  TestOutlierFinder tof{5_mm};
  kfOptions.extensions.outlierFinder
      .connect<&TestOutlierFinder::operator()<VectorMultiTrajectory>>(&tof);

  tester.test_ZeroFieldWithOutliers(kfZero, kfOptions, start, rng, false, false,
                                    true);
}

// TODO this is not really Kalman fitter specific. is probably better tested
// with a synthetic trajectory.
BOOST_AUTO_TEST_CASE(GlobalCovariance) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_GlobalCovariance(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
