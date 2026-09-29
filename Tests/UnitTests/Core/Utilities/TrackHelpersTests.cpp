// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <memory>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using enum TrackStateFlag;

namespace ActsTests {

namespace {

template <typename TrackContainer, typename FlagsPerState>
auto createTestTrack(TrackContainer& tc, const FlagsPerState& flagsPerState) {
  auto t = tc.makeTrack();

  for (const auto& flags : flagsPerState) {
    auto ts = t.appendTrackState();
    for (auto f : flags) {
      ts.typeFlags().setUnchecked(f);
    }
  }

  return t;
}

template <typename TrackContainer>
auto createTestTrackState(TrackContainer& tc) {
  auto t = tc.makeTrack();

  auto ts = t.appendTrackState();

  ts.allocateCalibrated(Vector2::Zero(), SquareMatrix2::Identity());
  ts.setProjectorSubspaceIndices(std::array{eBoundLoc0, eBoundLoc1});

  ts.predicted() = BoundVector::Zero();
  ts.predicted()[eBoundLoc0] = 1.;
  ts.predicted()[eBoundLoc1] = 1.;
  ts.predictedCovariance() = BoundMatrix::Identity() * 1.;

  ts.filtered() = BoundVector::Zero();
  ts.filtered()[eBoundLoc0] = 0.5;
  ts.filtered()[eBoundLoc1] = 0.5;
  ts.filteredCovariance() = BoundMatrix::Identity() * 0.5;

  ts.smoothed() = BoundVector::Zero();
  ts.smoothed()[eBoundLoc0] = 0.1;
  ts.smoothed()[eBoundLoc1] = 0.1;
  ts.smoothedCovariance() = BoundMatrix::Identity() * 0.1;

  return ts;
}

/// Which parameters a track state of the extrapolation test track carries.
enum class ParameterSlot { none, predicted, filtered, smoothed };

/// Build a track of measurement states on planes at z = 100, 200, 300, all
/// pointing along +z on the z axis, so the first state is the one closest to a
/// target plane at the origin.
template <typename TrackContainer>
auto createExtrapolationTestTrack(TrackContainer& tc,
                                  const std::vector<ParameterSlot>& slots) {
  auto t = tc.makeTrack();

  double z = 100_mm;
  for (ParameterSlot slot : slots) {
    TrackStatePropMask mask = TrackStatePropMask::None;
    if (slot == ParameterSlot::predicted) {
      mask = TrackStatePropMask::Predicted;
    } else if (slot == ParameterSlot::filtered) {
      mask = TrackStatePropMask::Filtered;
    } else if (slot == ParameterSlot::smoothed) {
      mask = TrackStatePropMask::Smoothed;
    }

    auto ts = t.appendTrackState(mask);
    ts.typeFlags().setUnchecked(HasMeasurement);
    ts.setReferenceSurface(Surface::makeShared<PlaneSurface>(
        Transform3{Translation3{Vector3{0, 0, z}}},
        std::make_shared<RectangleBounds>(1_m, 1_m)));

    // loc0 = loc1 = 0 and theta = 0 puts the state on the z axis pointing
    // along +z, so the path length to the target is just -z
    BoundVector params = BoundVector::Zero();
    params[eBoundQOverP] = 1 / 1_GeV;
    if (slot == ParameterSlot::predicted) {
      ts.predicted() = params;
      ts.predictedCovariance() = BoundMatrix::Identity();
    } else if (slot == ParameterSlot::filtered) {
      ts.filtered() = params;
      ts.filteredCovariance() = BoundMatrix::Identity();
    } else if (slot == ParameterSlot::smoothed) {
      ts.smoothed() = params;
      ts.smoothedCovariance() = BoundMatrix::Identity();
    }

    z += 100_mm;
  }

  return t;
}

/// Target plane at the origin, reached by the test track above.
std::shared_ptr<PlaneSurface> makeTargetSurface() {
  return Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(1_m, 1_m));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(CalculateQuantities) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto t = createTestTrack(tc, std::vector<std::vector<TrackStateFlag>>{
                                   {HasMeasurement},
                                   {HasMeasurement, IsOutlier},
                                   {HasMeasurement, IsSharedHit},
                                   {IsHole},
                                   {HasMeasurement, IsOutlier},
                                   {IsHole},
                                   {HasMeasurement, IsSharedHit},
                                   {HasMeasurement, IsOutlier},
                               });

  calculateTrackQuantities(t);

  BOOST_CHECK_EQUAL(t.nHoles(), 2);
  BOOST_CHECK_EQUAL(t.nMeasurements(), 3);
  BOOST_CHECK_EQUAL(t.nOutliers(), 3);
  BOOST_CHECK_EQUAL(t.nSharedHits(), 2);
}

BOOST_AUTO_TEST_CASE(TrimTrack) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto t = createTestTrack(tc, std::vector<std::vector<TrackStateFlag>>{
                                   {},
                                   {IsHole},
                                   {HasMeasurement},
                                   {HasMeasurement, IsOutlier},
                                   {HasMeasurement, IsSharedHit},
                                   {IsHole},
                                   {HasMeasurement, IsOutlier},
                                   {IsHole},
                                   {HasMeasurement},
                                   {HasMeasurement, IsOutlier},
                                   {},
                               });

  calculateTrackQuantities(t);

  BOOST_CHECK_EQUAL(t.nTrackStates(), 11);
  BOOST_CHECK_EQUAL(t.nHoles(), 3);
  BOOST_CHECK_EQUAL(t.nMeasurements(), 3);
  BOOST_CHECK_EQUAL(t.nOutliers(), 3);
  BOOST_CHECK_EQUAL(t.nSharedHits(), 1);

  trimTrackFront(t, true, true, true, true);
  calculateTrackQuantities(t);

  BOOST_CHECK_EQUAL(t.nTrackStates(), 9);
  BOOST_CHECK_EQUAL(t.nHoles(), 2);
  BOOST_CHECK_EQUAL(t.nMeasurements(), 3);
  BOOST_CHECK_EQUAL(t.nOutliers(), 3);
  BOOST_CHECK_EQUAL(t.nSharedHits(), 1);

  trimTrackBack(t, true, true, true, true);
  calculateTrackQuantities(t);

  BOOST_CHECK_EQUAL(t.nTrackStates(), 7);
  BOOST_CHECK_EQUAL(t.nHoles(), 2);
  BOOST_CHECK_EQUAL(t.nMeasurements(), 3);
  BOOST_CHECK_EQUAL(t.nOutliers(), 2);
  BOOST_CHECK_EQUAL(t.nSharedHits(), 1);

  trimTrack(t, true, true, true, true);
  calculateTrackQuantities(t);

  BOOST_CHECK_EQUAL(t.nTrackStates(), 7);
  BOOST_CHECK_EQUAL(t.nHoles(), 2);
  BOOST_CHECK_EQUAL(t.nMeasurements(), 3);
  BOOST_CHECK_EQUAL(t.nOutliers(), 2);
  BOOST_CHECK_EQUAL(t.nSharedHits(), 1);
}

BOOST_AUTO_TEST_CASE(CalculatePredictedChi2) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto ts = createTestTrackState(tc);

  // reference found by running the code
  BOOST_CHECK_CLOSE(calculatePredictedChi2(ts), 1., 1e-6);
}

BOOST_AUTO_TEST_CASE(CalculateFilteredChi2) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto ts = createTestTrackState(tc);

  // reference found by running the code
  BOOST_CHECK_CLOSE(calculateFilteredChi2(ts), 1., 1e-6);
}

BOOST_AUTO_TEST_CASE(CalculateSmoothedChi2) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto ts = createTestTrackState(tc);

  // reference found by running the code
  BOOST_CHECK_CLOSE(calculateSmoothedChi2(ts), 1. / 45., 1e-6);
}

BOOST_AUTO_TEST_CASE(CalculateUnbiasedParametersCovariance) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto ts = createTestTrackState(tc);

  auto [params, cov] =
      calculateUnbiasedParametersCovariance(AnyConstTrackStateProxy{ts});

  // reference found by running the code
  BoundVector refParams = BoundVector::Zero();
  refParams[eBoundLoc0] = 1. / 9.;
  refParams[eBoundLoc1] = 1. / 9.;
  BoundMatrix refCov = BoundMatrix::Identity() * 0.1;
  refCov(eBoundLoc0, eBoundLoc0) = 1. / 9.;
  refCov(eBoundLoc1, eBoundLoc1) = 1. / 9.;

  CHECK_CLOSE_ABS(params, refParams, 1e-6);
  CHECK_CLOSE_ABS(cov, refCov, 1e-6);
}

BOOST_AUTO_TEST_CASE(FindTrackStateForExtrapolationStrategies) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto target = makeTargetSurface();
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  auto t = createExtrapolationTestTrack(
      tc, {ParameterSlot::smoothed, ParameterSlot::smoothed,
           ParameterSlot::smoothed});

  auto first = findTrackStateForExtrapolation(
      gctx, t, *target, TrackExtrapolationStrategy::first);
  BOOST_REQUIRE(first.ok());
  CHECK_CLOSE_ABS(first->second, -100_mm, 1e-6);

  auto last = findTrackStateForExtrapolation(gctx, t, *target,
                                             TrackExtrapolationStrategy::last);
  BOOST_REQUIRE(last.ok());
  CHECK_CLOSE_ABS(last->second, -300_mm, 1e-6);

  // the first state is the closer one, so `firstOrLast` has to agree with
  // `first`
  auto firstOrLast = findTrackStateForExtrapolation(
      gctx, t, *target, TrackExtrapolationStrategy::firstOrLast);
  BOOST_REQUIRE(firstOrLast.ok());
  BOOST_CHECK_EQUAL(firstOrLast->first.index(), first->first.index());
  CHECK_CLOSE_ABS(firstOrLast->second, -100_mm, 1e-6);
}

BOOST_AUTO_TEST_CASE(FindTrackStateForExtrapolationParameterSlots) {
  auto target = makeTargetSurface();
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  // a state carrying only predicted is usable, since that is what
  // `TrackProxy::createParametersFromState` would start from
  for (ParameterSlot slot : {ParameterSlot::predicted, ParameterSlot::filtered,
                             ParameterSlot::smoothed}) {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto t = createExtrapolationTestTrack(tc, {slot, slot, slot});

    auto result = findTrackStateForExtrapolation(
        gctx, t, *target, TrackExtrapolationStrategy::firstOrLast);
    BOOST_REQUIRE(result.ok());
    CHECK_CLOSE_ABS(result->second, -100_mm, 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(FindTrackStateForExtrapolationPartialParameters) {
  auto target = makeTargetSurface();
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  // a seed track only carries parameters on the bottom space point, so
  // `firstOrLast` cannot intersect the last state and has to fall back to the
  // first one
  {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto t = createExtrapolationTestTrack(
        tc,
        {ParameterSlot::predicted, ParameterSlot::none, ParameterSlot::none});

    auto result = findTrackStateForExtrapolation(
        gctx, t, *target, TrackExtrapolationStrategy::firstOrLast);
    BOOST_REQUIRE(result.ok());
    CHECK_CLOSE_ABS(result->second, -100_mm, 1e-6);
  }

  // and the other way around, where the fallback has to pick the state that is
  // further away
  {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto t = createExtrapolationTestTrack(
        tc,
        {ParameterSlot::none, ParameterSlot::none, ParameterSlot::filtered});

    auto result = findTrackStateForExtrapolation(
        gctx, t, *target, TrackExtrapolationStrategy::firstOrLast);
    BOOST_REQUIRE(result.ok());
    CHECK_CLOSE_ABS(result->second, -300_mm, 1e-6);
  }

  // asking for the end that carries no parameters is an error rather than a
  // read of an unallocated parameter slot
  {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto t = createExtrapolationTestTrack(
        tc,
        {ParameterSlot::filtered, ParameterSlot::none, ParameterSlot::none});

    auto result = findTrackStateForExtrapolation(
        gctx, t, *target, TrackExtrapolationStrategy::last);
    BOOST_REQUIRE(!result.ok());
    BOOST_CHECK_EQUAL(result.error(),
                      TrackExtrapolationError::CompatibleTrackStateNotFound);
  }

  // no end carries parameters
  {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto t = createExtrapolationTestTrack(
        tc, {ParameterSlot::none, ParameterSlot::none, ParameterSlot::none});

    auto result = findTrackStateForExtrapolation(
        gctx, t, *target, TrackExtrapolationStrategy::firstOrLast);
    BOOST_REQUIRE(!result.ok());
    BOOST_CHECK_EQUAL(result.error(),
                      TrackExtrapolationError::CompatibleTrackStateNotFound);
  }

  // a track without measurement states at all
  {
    TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto t = tc.makeTrack();

    auto result = findTrackStateForExtrapolation(
        gctx, t, *target, TrackExtrapolationStrategy::firstOrLast);
    BOOST_REQUIRE(!result.ok());
    BOOST_CHECK_EQUAL(result.error(),
                      TrackExtrapolationError::CompatibleTrackStateNotFound);
  }
}

BOOST_AUTO_TEST_CASE(FindTrackStateForExtrapolationUnreachable) {
  TrackContainer tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  auto t = createExtrapolationTestTrack(
      tc, {ParameterSlot::smoothed, ParameterSlot::smoothed,
           ParameterSlot::smoothed});

  // the track runs along +z on the z axis, so a target plane offset in x is
  // outside the bounds it can reach
  auto target = Surface::makeShared<PlaneSurface>(
      Transform3{Translation3{Vector3{10_m, 0, 0}}},
      std::make_shared<RectangleBounds>(1_mm, 1_mm));

  auto result = findTrackStateForExtrapolation(
      gctx, t, *target, TrackExtrapolationStrategy::firstOrLast);
  BOOST_REQUIRE(!result.ok());
  BOOST_CHECK_EQUAL(result.error(),
                    TrackExtrapolationError::ReferenceSurfaceUnreachable);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
