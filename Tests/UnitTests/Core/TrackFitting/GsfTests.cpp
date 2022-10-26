// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GaussianSumFitter.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"

#include <algorithm>
#include <memory>
#include <random>

#include "FitterTestsCommon.hpp"

namespace {

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::UnitLiterals;
using namespace Acts::Experimental;

Acts::GainMatrixUpdater kfUpdater;

GsfExtensions<VectorMultiTrajectory> getExtensions() {
  GsfExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  extensions.updater
      .connect<&Acts::GainMatrixUpdater::operator()<VectorMultiTrajectory>>(
          &kfUpdater);
  return extensions;
}

FitterTester tester;

const auto logger = getDefaultLogger("GSF", Logging::INFO);

using Stepper = Acts::MultiEigenStepperLoop<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;

auto gsfZeroPropagator =
    makeConstantFieldPropagator<Stepper>(tester.geometry, 0_T);
const GaussianSumFitter<Propagator, VectorMultiTrajectory> gsfZero(
    std::move(gsfZeroPropagator));

std::default_random_engine rng(42);

auto makeDefaultGsfOptions() {
  return GsfOptions<VectorMultiTrajectory>{
      tester.geoCtx,   tester.magCtx,          tester.calCtx,
      getExtensions(), LoggerWrapper{*logger}, PropagatorPlainOptions()};
}

// A Helper type to allow us to put the MultiComponentBoundTrackParameters into
// the function so that it can also be used as SingleBoundTrackParameters for
// the MeasurementsCreator
template <typename charge_t>
struct MultiCmpsParsInterface : public SingleBoundTrackParameters<charge_t> {
  MultiComponentBoundTrackParameters<charge_t> multi_pars;

  MultiCmpsParsInterface(const MultiComponentBoundTrackParameters<charge_t> &p)
      : SingleBoundTrackParameters<charge_t>(
            p.referenceSurface().getSharedPtr(), p.parameters(),
            p.covariance()),
        multi_pars(p) {}

  operator MultiComponentBoundTrackParameters<charge_t>() const {
    return multi_pars;
  }
};

auto makeParameters() {
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
  Acts::CurvilinearTrackParameters cp(mPos4, 0_degree, 90_degree, 1_GeV, 1_e,
                                      cov);

  // Construct bound multi component parameters from curvilinear ones
  Acts::BoundVector deltaLOC0 = Acts::BoundVector::Zero();
  deltaLOC0[eBoundLoc0] = 0.5_mm;

  Acts::BoundVector deltaLOC1 = Acts::BoundVector::Zero();
  deltaLOC1[eBoundLoc1] = 0.5_mm;

  Acts::BoundVector deltaQOP = Acts::BoundVector::Zero();
  deltaQOP[eBoundQOverP] = 0.01_GeV;

  std::vector<std::tuple<double, BoundVector, BoundSymMatrix>> cmps = {
      {0.2, cp.parameters(), cov},
      {0.2, cp.parameters() + deltaLOC0 + deltaLOC1 + deltaQOP, cov},
      {0.2, cp.parameters() + deltaLOC0 - deltaLOC1 - deltaQOP, cov},
      {0.2, cp.parameters() - deltaLOC0 + deltaLOC1 + deltaQOP, cov},
      {0.2, cp.parameters() - deltaLOC0 - deltaLOC1 - deltaQOP, cov}};

  return MultiCmpsParsInterface<SinglyCharged>(
      Acts::MultiComponentBoundTrackParameters<SinglyCharged>(
          cp.referenceSurface().getSharedPtr(), cmps));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFittingGsf)

BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldNoSurfaceForward(gsfZero, options, multi_pars, rng, true,
                                        true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceForward) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldWithSurfaceForward(gsfZero, options, multi_pars, rng,
                                          true, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceBackward) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldWithSurfaceBackward(gsfZero, options, multi_pars, rng,
                                           true, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceAtExit) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldWithSurfaceBackward(gsfZero, options, multi_pars, rng,
                                           true, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldShuffled) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldShuffled(gsfZero, options, multi_pars, rng, true, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithHole) {
  auto options = makeDefaultGsfOptions();
  auto multi_pars = makeParameters();

  tester.test_ZeroFieldWithHole(gsfZero, options, multi_pars, rng, true, true);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithOutliers) {
  // fitter options w/o target surface. outlier distance is set to be below the
  // default outlier distance in the `MeasurementsCreator`
  TestOutlierFinder tof{5_mm};
  auto options = makeDefaultGsfOptions();
  options.extensions.outlierFinder
      .connect<&TestOutlierFinder::operator()<VectorMultiTrajectory>>(&tof);

  auto multi_pars = makeParameters();

  tester.test_ZeroFieldWithOutliers(gsfZero, options, multi_pars, rng, true,
                                    true);
}

// NOTE This test makes no sense for the GSF since there is always reverse
// filtering BOOST_AUTO_TEST_CASE(ZeroFieldWithReverseFiltering) { ... }

// TODO this is not really Kalman fitter specific. is probably better tested
// with a synthetic trajectory.
BOOST_AUTO_TEST_CASE(GlobalCovariance) {
  auto options = makeDefaultGsfOptions();
  auto multi_pars = makeParameters();

  tester.test_GlobalCovariance(gsfZero, options, multi_pars, rng);
}

BOOST_AUTO_TEST_SUITE_END()
