// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/GaussianSumFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"

#include <algorithm>
#include <memory>
#include <random>

#include "FitterTestsCommon.hpp"

namespace {

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::UnitLiterals;

using Stepper = Acts::MultiEigenStepperLoop<>;
using Propagator =
    Acts::Propagator<Stepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using Gsf = Acts::KalmanFitter<Propagator>;

KalmanUpdater kfUpdater;
KalmanSmoother kfSmoother;

KalmanFitterExtensions getExtensions() {
  KalmanFitterExtensions extensions;
  extensions.calibrator.connect<&testSourceLinkCalibrator>();
  extensions.updater.connect<&KalmanUpdater::operator()>(&kfUpdater);
  extensions.smoother.connect<&KalmanSmoother::operator()>(&kfSmoother);
  return extensions;
}

// reconstruction propagator and fitter
const auto logger = getDefaultLogger("GSF", Logging::INFO);
const auto gsfZeroPropagator =
    makeConstantFieldPropagator<Stepper>(geometry, 0_T);
const auto gsfZero = GaussianSumFitter(std::move(gsfZeroPropagator));

std::default_random_engine rng(42);

auto makeDefaultGsfOptions() {
  return GsfOptions{geoCtx, magCtx, calCtx, getExtensions(),
                             LoggerWrapper{*logger},
                             PropagatorPlainOptions()};
}

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFittingKalmanFitter)

BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldNoSurfaceForward(gsfZero, options, rng);
}

/*
BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceForward) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldWithSurfaceForward(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceBackward) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldWithSurfaceBackward(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceAtExit) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldWithSurfaceAtExit(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldShuffled) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldShuffled(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithHole) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldWithHole(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithOutliers) {
  // fitter options w/o target surface. outlier distance is set to be below the
  // default outlier distance in the `MeasurementsCreator`
  auto extensions = getExtensions();
  TestOutlierFinder tof{5_mm};
  extensions.outlierFinder.connect<&TestOutlierFinder::operator()>(&tof);

  KalmanFitterOptions options(geoCtx, magCtx, calCtx, extensions,
                                LoggerWrapper{*kfLogger},
                                PropagatorPlainOptions());

  test_ZeroFieldWithOutliers(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithReverseFiltering) {
  auto test = [](double threshold, bool reverse, bool expected_reversed,
                 bool expected_smoothed) {
    // Reverse filtering threshold set at 0.5 GeV
    auto extensions = getExtensions();
    TestReverseFilteringLogic trfl{threshold};
    extensions.reverseFilteringLogic
        .connect<&TestReverseFilteringLogic::operator()>(&trfl);

    KalmanFitterOptions options(geoCtx, magCtx, calCtx, extensions,
                                  LoggerWrapper{*kfLogger},
                                  PropagatorPlainOptions());
    options.reversedFiltering = reverse;
    test_ZeroFieldWithReverseFiltering(kfZero, options, rng,
                                       expected_reversed, expected_smoothed);
  };

  // Track of 1 GeV with a threshold set at 0.1 GeV, reversed filtering should
  // not be used
  test(0.1_GeV, false, false, true);

  // Track of 1 GeV with a threshold set at 10 GeV, reversed filtering should
  // be used
  test(10._GeV, false, true, false);

  // Track of 1 GeV with a threshold set at 10 GeV, reversed filtering should
  // be used
  test(0.1_GeV, true, true, false);
}

// TODO this is not really Kalman fitter specific. is probably better tested
// with a synthetic trajectory.
BOOST_AUTO_TEST_CASE(GlobalCovariance) {
  auto options = makeDefaultGsfOptions();

  test_GlobalCovariance(kfZero, options, rng);
}
*/

BOOST_AUTO_TEST_SUITE_END()

