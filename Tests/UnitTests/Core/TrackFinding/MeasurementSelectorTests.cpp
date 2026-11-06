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
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(MeasurementSelectorTests)

namespace {

const auto gctx = GeometryContext::dangerouslyDefaultConstruct();

std::shared_ptr<Surface> makeSurface() {
  auto bounds = std::make_shared<RectangleBounds>(100., 100.);
  auto surface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), bounds);
  surface->assignGeometryId(
      GeometryIdentifier().withVolume(1).withLayer(1).withSensitive(1));
  return surface;
}

/// Adds a track state candidate on @p surface.
///
/// A fit that diverges ends up with non-finite predicted parameters, and
/// inverting a sufficiently inflated predicted covariance can overflow and
/// yield a non-finite chi2 as well. The non-finite value is injected directly
/// rather than provoked through an overflow, so that the chi2 is non-finite on
/// every platform instead of depending on how the matrix inverse rounds.
VectorMultiTrajectory::TrackStateProxy addCandidate(
    VectorMultiTrajectory& traj, const std::shared_ptr<Surface>& surface,
    bool diverged) {
  auto ts = traj.makeTrackState(TrackStatePropMask::All);
  ts.setReferenceSurface(surface);

  ts.predicted() = BoundVector::Zero();
  // Avoid theta = 0, whose pseudorapidity is infinite, since the selector bins
  // its cuts in eta.
  ts.predicted()[eBoundTheta] = std::numbers::pi / 2.;
  if (diverged) {
    ts.predicted()[eBoundLoc0] = std::numeric_limits<double>::quiet_NaN();
  }
  ts.predictedCovariance() = BoundMatrix::Identity() * 1e-2;

  ts.allocateCalibrated(2);
  ts.calibrated<2>() << 0.0, 0.0;
  ts.calibratedCovariance<2>() = SquareMatrix<2>::Identity() * 1e-2;
  ts.setProjectorSubspaceIndices(
      std::array<std::uint8_t, 2>{eBoundLoc0, eBoundLoc1});
  return ts;
}

}  // namespace

/// Every comparison against NaN is false, so a candidate with a non-finite
/// chi2 used to leave minIndex unset while still being counted as passing the
/// chi2 cut. The selector then returned `candidates.begin() - 1`, which the
/// caller dereferences.
BOOST_AUTO_TEST_CASE(NonFiniteChi2IsRejected) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  auto surface = makeSurface();

  VectorMultiTrajectory traj;
  std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{
      addCandidate(traj, surface, true)};

  MeasurementSelectorCuts cuts;
  cuts.chi2CutOff = {15};
  cuts.numMeasurementsCutOff = {1};
  MeasurementSelector selector{cuts};

  bool isOutlier = false;
  auto result = selector.select<VectorMultiTrajectory>(gctx, candidates,
                                                       isOutlier, *logger);

  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE(!std::isfinite(candidates.front().chi2()));

  auto [begin, end] = *result;
  // The returned range must be a valid subrange of the candidates.
  BOOST_CHECK_GE(begin - candidates.begin(), 0);
  BOOST_CHECK_LE(end - candidates.begin(),
                 static_cast<std::ptrdiff_t>(candidates.size()));
  // The only candidate is unusable, so nothing may be selected.
  BOOST_CHECK_EQUAL(end - begin, 0);
}

/// A candidate with a finite chi2 must still be selected when a non-finite one
/// is present alongside it, rather than the non-finite one being preferred.
BOOST_AUTO_TEST_CASE(FiniteChi2StillSelectedAlongsideNonFinite) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  auto surface = makeSurface();

  VectorMultiTrajectory traj;
  std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{
      addCandidate(traj, surface, true), addCandidate(traj, surface, false)};

  MeasurementSelectorCuts cuts;
  cuts.chi2CutOff = {15};
  cuts.numMeasurementsCutOff = {1};
  MeasurementSelector selector{cuts};

  bool isOutlier = false;
  auto result = selector.select<VectorMultiTrajectory>(gctx, candidates,
                                                       isOutlier, *logger);

  BOOST_REQUIRE(result.ok());
  auto [begin, end] = *result;
  BOOST_CHECK_GE(begin - candidates.begin(), 0);
  BOOST_REQUIRE_EQUAL(end - begin, 1);
  BOOST_CHECK(std::isfinite(begin->chi2()));
}

BOOST_AUTO_TEST_CASE(PhiSeamResidualAndThetaDependentCuts) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  const auto surface = makeSurface();
  MeasurementSelectorCuts cuts;
  // Preserve the existing direct theta comparison of these bin edges.
  cuts.etaBins = {1.};
  cuts.chi2CutOff = {.1, 1.};
  cuts.chi2CutOffOutlier = {0.};
  cuts.numMeasurementsCutOff = {1, 1};
  MeasurementSelector selector(cuts);

  for (double theta : {std::numbers::pi / 2., .3}) {
    VectorMultiTrajectory trajectory;
    auto candidate = addCandidate(trajectory, surface, false);
    candidate.predicted()[eBoundDir0] = std::numbers::pi - .01;
    candidate.predicted()[eBoundDir1] = theta;
    candidate.setProjectorSubspaceIndices(
        std::array<std::uint8_t, 2>{eBoundDir0, eBoundLoc0});
    candidate.calibrated<2>() << -std::numbers::pi + .01, .1;
    std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{candidate};
    bool isOutlier = false;
    const auto result = selector.select<VectorMultiTrajectory>(
        gctx, candidates, isOutlier, *logger);
    BOOST_REQUIRE(result.ok());
    // Variance .02 in each measured slot; residuals (.02, .1).
    BOOST_CHECK_SMALL(candidate.chi2() - .52, 1e-6);
    BOOST_CHECK_EQUAL(result->second - result->first, theta == .3 ? 0 : 1);
    BOOST_CHECK(!isOutlier);
  }
}

BOOST_AUTO_TEST_CASE(PhiSeamOutlierAndNonfiniteDirection) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  const auto surface = makeSurface();
  VectorMultiTrajectory trajectory;
  auto invalid = addCandidate(trajectory, surface, false);
  invalid.predicted()[eBoundDir0] = std::numeric_limits<double>::quiet_NaN();
  auto valid = addCandidate(trajectory, surface, false);
  valid.predicted()[eBoundDir0] = std::numbers::pi - .1;
  valid.setProjectorSubspaceIndices(
      std::array<std::uint8_t, 2>{eBoundDir0, eBoundLoc0});
  valid.calibrated<2>() << -std::numbers::pi + .1, .1;
  std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{invalid,
                                                                 valid};
  MeasurementSelectorCuts cuts;
  cuts.chi2CutOff = {1.};
  cuts.chi2CutOffOutlier = {3.};
  MeasurementSelector selector(cuts);
  bool isOutlier = false;
  const auto result = selector.select<VectorMultiTrajectory>(
      gctx, candidates, isOutlier, *logger);
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->second - result->first, 1);
  BOOST_CHECK_EQUAL(result->first->index(), valid.index());
  BOOST_CHECK_SMALL(valid.chi2() - 2.5, 1e-6);
  BOOST_CHECK(!std::isfinite(invalid.chi2()));
  BOOST_CHECK(isOutlier);

  // Position residuals must not wrap like direction coordinates.
  valid.calibrated<2>()[1] = 2. * std::numbers::pi + .1;
  const auto rejected = selector.select<VectorMultiTrajectory>(
      gctx, candidates, isOutlier, *logger);
  BOOST_REQUIRE(rejected.ok());
  BOOST_CHECK(rejected->first == rejected->second);
  BOOST_CHECK(!isOutlier);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
