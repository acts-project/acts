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
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"
#include "Acts/TrackFitting/ReferenceTrajectoryBuilder.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <random>
#include <ranges>

#include "FitterTestsCommon.hpp"

using namespace Acts;
using namespace Acts::detail::Test;
using namespace Acts::UnitLiterals;

namespace ActsTests {

using ConstantFieldStepper = EigenStepper<>;
using ConstantFieldPropagator = Propagator<ConstantFieldStepper, Navigator>;
using ReferenceTrajectoryBuilder =
    Experimental::ReferenceTrajectoryBuilder<ConstantFieldPropagator,
                                             VectorMultiTrajectory>;

// Construct initial track parameters.
BoundTrackParameters makeParameters(
    const ParticleHypothesis& particleHypothesis) {
  // create covariance matrix from reasonable standard deviations
  BoundVector stddev;
  stddev[eBoundLoc0] = 100_um;
  stddev[eBoundLoc1] = 100_um;
  stddev[eBoundTime] = 25_ns;
  stddev[eBoundPhi] = 2_degree;
  stddev[eBoundTheta] = 2_degree;
  stddev[eBoundQOverP] = 1 / 100_GeV;
  BoundMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  Vector4 mPos4(-3_m, 0., 0., 42_ns);
  return BoundTrackParameters::createCurvilinear(
      mPos4, 0_degree, 90_degree, 1_e / 1_GeV, cov, particleHypothesis);
}

// Instantiate the tester
const FitterTester tester;

/// Adapter exposing the reference trajectory builder through the `fit`
/// interface used by FitterTester: build, attach and calibrate measurements,
/// filter and smooth.
struct ReferenceTrajectoryFitter {
  ReferenceTrajectoryBuilder builder;
  GainMatrixUpdater updater;
  MbfSmoother smoother;

  template <typename source_link_iterator_t, typename track_container_t>
  Result<typename track_container_t::TrackProxy> fit(
      source_link_iterator_t begin, source_link_iterator_t end,
      const BoundTrackParameters& start,
      const ReferenceTrajectoryBuilder::Options& options,
      track_container_t& tracks) const {
    auto buildResult = builder.build(start, options, tracks);
    if (!buildResult.ok()) {
      return buildResult.error();
    }
    auto track = *buildResult;

    SourceLinkSurfaceAccessor surfaceAccessor;
    surfaceAccessor.connect<&TestSourceLink::SurfaceAccessor::operator()>(
        &tester.surfaceAccessor);
    builder.attachSourceLinks(track, std::ranges::subrange(begin, end),
                              surfaceAccessor);

    ReferenceTrajectoryBuilder::Calibrator calibrator;
    calibrator.connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
    builder.calibrateMeasurements(options.geoContext, tester.calCtx, track,
                                  calibrator);

    ReferenceTrajectoryBuilder::Updater updaterDelegate;
    updaterDelegate
        .connect<&GainMatrixUpdater::operator()<VectorMultiTrajectory>>(
            &updater);
    if (auto res = builder.filter(options.geoContext, track, updaterDelegate);
        !res.ok()) {
      return res.error();
    }

    if (auto res = smoother(options.geoContext, tracks.trackStateContainer(),
                            track.tipIndex());
        !res.ok()) {
      return res.error();
    }

    calculateTrackQuantities(track);

    return track;
  }
};

const ReferenceTrajectoryFitter rtFitter{
    ReferenceTrajectoryBuilder(
        makeConstantFieldPropagator<ConstantFieldStepper>(tester.geometry,
                                                          0_T)),
    {},
    {}};

std::default_random_engine rng(42);

ReferenceTrajectoryBuilder::Options makeDefaultOptions() {
  return ReferenceTrajectoryBuilder::Options(
      tester.geoCtx, tester.magCtx,
      PropagatorPlainOptions(tester.geoCtx, tester.magCtx));
}

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

BOOST_AUTO_TEST_CASE(OutputTrackParticleHypothesis) {
  auto start = makeParameters(ParticleHypothesis::muon());
  auto options = makeDefaultOptions();

  tester.test_ParticleHypothesis(rtFitter, options, start, rng);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
