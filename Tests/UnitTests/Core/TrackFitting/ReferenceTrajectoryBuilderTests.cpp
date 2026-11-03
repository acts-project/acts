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
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"
#include "Acts/TrackFitting/ReferenceTrajectoryBuilder.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <numbers>
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
  stddev[eBoundDir0] = 2_degree;
  stddev[eBoundDir1] = 2_degree;
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

namespace {
using Builder = Experimental::ReferenceTrajectoryBuilder<
    Propagator<StraightLineStepper, DirectNavigator>, VectorMultiTrajectory>;

struct CorrectionUpdater {
  BoundVector delta = BoundVector::Zero();
  Result<void> update(const GeometryContext& gctx,
                      VectorMultiTrajectory::TrackStateProxy state,
                      const Logger&) const {
    const auto result =
        state.referenceSurface().boundParametrization(gctx).canonicalize(
            state.predicted() + delta, state.predictedCovariance());
    if (!result.ok()) {
      return result.error();
    }
    state.filtered() = result->first;
    state.filteredCovariance() = result->second;
    return Result<void>::success();
  }
};

struct Fixture {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  TrackContainer<VectorTrackContainer, VectorMultiTrajectory, ValueHolder>
      tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  Builder builder{Propagator{StraightLineStepper{}, DirectNavigator{}}};
  std::shared_ptr<PlaneSurface> surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity());

  template <typename Track>
  auto append(Track& track, double phi, double theta) {
    auto state = track.appendTrackState(TrackStatePropMask::All);
    state.setReferenceSurface(surface);
    state.predicted() << 0., 0., phi, theta, .1, 0.;
    state.predictedCovariance().setIdentity();
    state.predictedCovariance()(eBoundLoc0, eBoundDir1) = .25;
    state.predictedCovariance()(eBoundDir1, eBoundLoc0) = .25;
    state.jacobian().setIdentity();
    return state;
  }
};
}  // namespace

BOOST_AUTO_TEST_SUITE(ReferenceTrajectoryBuilderSuite)

BOOST_FIXTURE_TEST_CASE(ReferenceResidualAcrossPhiSeam, Fixture) {
  auto track = tracks.makeTrack();
  auto first = append(track, std::numbers::pi - .1, 1.);
  first.typeFlags().setHasMeasurement();
  auto second = append(track, std::numbers::pi - .05, 1.);
  auto third = append(track, 0., 1.);
  third.jacobian()(eBoundLoc0, eBoundDir0) = 2.;
  track.linkForward();

  CorrectionUpdater updater;
  updater.delta[eBoundDir0] = .2;
  Builder::Updater delegate;
  delegate.connect<&CorrectionUpdater::update>(&updater);
  BOOST_REQUIRE(builder.filter(gctx, track, delegate).ok());
  CHECK_CLOSE_ABS(second.predicted()[eBoundDir0], -std::numbers::pi + .15,
                  1e-12);
  CHECK_CLOSE_ABS(third.predicted()[eBoundDir0], .2, 1e-12);
  CHECK_CLOSE_ABS(third.predicted()[eBoundLoc0], .4, 1e-12);
  const BoundMatrix expected = third.jacobian() * first.filteredCovariance() *
                               third.jacobian().transpose();
  CHECK_CLOSE_ABS(third.predictedCovariance(), expected, 1e-12);
}

BOOST_FIXTURE_TEST_CASE(PredictionFoldTransformsCovarianceAndJacobian,
                        Fixture) {
  auto track = tracks.makeTrack();
  auto first = append(track, .2, .5);
  first.typeFlags().setHasMeasurement();
  auto second = append(track, .2, .1);
  track.linkForward();

  CorrectionUpdater updater;
  updater.delta[eBoundDir1] = -.2;
  Builder::Updater delegate;
  delegate.connect<&CorrectionUpdater::update>(&updater);
  BOOST_REQUIRE(builder.filter(gctx, track, delegate).ok());
  CHECK_CLOSE_ABS(second.predicted()[eBoundDir1], .1, 1e-12);
  CHECK_CLOSE_ABS(second.predictedCovariance()(eBoundLoc0, eBoundDir1), -.25,
                  1e-12);
  CHECK_CLOSE_ABS(second.jacobian()(eBoundDir1, eBoundDir1), -1., 1e-12);
  CHECK_CLOSE_ABS(second.filtered(), second.predicted(), 1e-12);
}

BOOST_FIXTURE_TEST_CASE(PoleFailurePreservesPredictionAndStops, Fixture) {
  auto track = tracks.makeTrack();
  auto first = append(track, .2, .5);
  first.typeFlags().setHasMeasurement();
  auto second = append(track, .2, .25);
  auto third = append(track, .2, .5);
  const BoundVector reference = second.predicted();
  const BoundMatrix covariance = second.predictedCovariance();
  track.linkForward();

  CorrectionUpdater updater;
  updater.delta[eBoundDir1] = -.25;
  Builder::Updater delegate;
  delegate.connect<&CorrectionUpdater::update>(&updater);
  const auto result = builder.filter(gctx, track, delegate);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK(result.error() ==
              BoundParametrizationError::UndefinedCanonicalization);
  CHECK_CLOSE_ABS(second.predicted(), reference, 1e-12);
  CHECK_CLOSE_ABS(second.predictedCovariance(), covariance, 1e-12);
  CHECK_CLOSE_ABS(second.jacobian(), BoundMatrix::Identity(), 1e-12);
  CHECK_CLOSE_ABS(third.predicted()[eBoundDir1], .5, 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
