// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/TrackFitting/ReferenceTrajectoryBuilder.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <numbers>

using namespace Acts;

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
