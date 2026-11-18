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
#include "Acts/EventData/BoundParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <stdexcept>
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

/// Adds a track state candidate on @p surface with a measurement of dimension
/// @p dim and the chi2 @p chi2 with respect to the prediction.
template <std::size_t dim>
VectorMultiTrajectory::TrackStateProxy addCandidateWithChi2(
    VectorMultiTrajectory& traj, const std::shared_ptr<Surface>& surface,
    double chi2) {
  auto ts = traj.makeTrackState(TrackStatePropMask::All);
  ts.setReferenceSurface(surface);

  ts.predicted() = BoundVector::Zero();
  ts.predicted()[eBoundDir1] = std::numbers::pi / 2.;
  // The predicted covariance in the measurement subspace is zero, so the chi2
  // is the squared residual
  ts.predictedCovariance() = BoundMatrix::Zero();

  ts.allocateCalibrated(dim);
  ts.template calibrated<dim>() = Vector<dim>::Zero();
  ts.template calibrated<dim>()[0] = std::sqrt(chi2);
  ts.template calibratedCovariance<dim>() = SquareMatrix<dim>::Identity();
  std::array<std::uint8_t, dim> indices{};
  for (std::size_t i = 0; i < dim; ++i) {
    indices[i] = static_cast<std::uint8_t>(i);
  }
  ts.setProjectorSubspaceIndices(indices);
  return ts;
}

/// Outcome of the selection of a single candidate
enum class Selection { Measurement, Outlier, Rejected };

/// Run the selector on a single candidate and return the outcome.
template <std::size_t dim>
Selection selectSingle(const MeasurementSelector& selector, double chi2) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  auto surface = makeSurface();

  VectorMultiTrajectory traj;
  std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{
      addCandidateWithChi2<dim>(traj, surface, chi2)};

  bool isOutlier = false;
  auto result = selector.select<VectorMultiTrajectory>(
      surface->boundParametrization(gctx), candidates, isOutlier, *logger);
  BOOST_REQUIRE(result.ok());
  auto [begin, end] = *result;
  BOOST_REQUIRE_LE(end - begin, 1);
  if (begin == end) {
    BOOST_CHECK(!isOutlier);
    return Selection::Rejected;
  }
  return isOutlier ? Selection::Outlier : Selection::Measurement;
}

/// Sets a 2D measurement at (@p loc0, 0) with variance 1e-2 on @p ts.
///
/// Together with the predicted covariance of `addCandidate`, a finite
/// measurement has the chi2 `loc0^2 / 0.02`.
void setMeasurement(VectorMultiTrajectory::TrackStateProxy& ts, double loc0) {
  ts.allocateCalibrated(2);
  ts.calibrated<2>() << loc0, 0.0;
  ts.calibratedCovariance<2>() = SquareMatrix<2>::Identity() * 1e-2;
  ts.setProjectorSubspaceIndices(
      std::array<std::uint8_t, 2>{eBoundLoc0, eBoundLoc1});
}

/// Adds a candidate sharing the prediction of @p first, like the track state
/// creator does for all but the first candidate on a surface.
VectorMultiTrajectory::TrackStateProxy addSharedCandidate(
    VectorMultiTrajectory& traj, VectorMultiTrajectory::TrackStateProxy& first,
    double loc0) {
  auto ts = traj.makeTrackState(TrackStatePropMask::Calibrated);
  ts.shareFrom(first, TrackStatePropMask::Predicted);
  ts.setReferenceSurface(first.referenceSurface().getSharedPtr());
  setMeasurement(ts, loc0);
  return ts;
}

MeasurementSelector makeSelector(std::size_t numMeasurements) {
  MeasurementSelectorCuts cuts;
  cuts.chi2CutOff = {15};
  cuts.numMeasurementsCutOff = {numMeasurements};
  return MeasurementSelector{cuts};
}

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // namespace

/// Without the per dimension flag every dimension uses the same chi2 cut.
BOOST_AUTO_TEST_CASE(Chi2CutSameForAllDimensions) {
  MeasurementSelectorCuts cuts;
  cuts.chi2CutOff = {15};
  cuts.chi2CutOffOutlier = {25};
  MeasurementSelector selector{cuts};

  BOOST_CHECK(selectSingle<1>(selector, 14.9) == Selection::Measurement);
  BOOST_CHECK(selectSingle<1>(selector, 15.1) == Selection::Outlier);
  BOOST_CHECK(selectSingle<2>(selector, 15.1) == Selection::Outlier);
  BOOST_CHECK(selectSingle<3>(selector, 15.1) == Selection::Outlier);
  BOOST_CHECK(selectSingle<3>(selector, 25.1) == Selection::Rejected);
}

/// With the per dimension flag the cuts are the 1D cuts for 1D measurements and
/// the chi2 values with the same p-value for higher dimensions.
BOOST_AUTO_TEST_CASE(Chi2CutPerDimension) {
  MeasurementSelectorCuts cuts;
  cuts.chi2CutOff = {15};
  cuts.chi2CutOffOutlier = {25};
  cuts.chi2CutPerDimension = true;
  MeasurementSelector selector{cuts};

  // Same p-value as chi2 = 15 (25) for one degree of freedom
  constexpr double cut2D = 18.28;
  constexpr double cut3D = 20.96;
  constexpr double outlier3D = 31.81;

  BOOST_CHECK(selectSingle<1>(selector, 14.9) == Selection::Measurement);
  BOOST_CHECK(selectSingle<1>(selector, 15.1) == Selection::Outlier);
  BOOST_CHECK(selectSingle<1>(selector, 25.1) == Selection::Rejected);

  BOOST_CHECK(selectSingle<2>(selector, cut2D - 0.1) == Selection::Measurement);
  BOOST_CHECK(selectSingle<2>(selector, cut2D + 0.1) == Selection::Outlier);

  BOOST_CHECK(selectSingle<3>(selector, cut3D - 0.1) == Selection::Measurement);
  BOOST_CHECK(selectSingle<3>(selector, cut3D + 0.1) == Selection::Outlier);
  BOOST_CHECK(selectSingle<3>(selector, outlier3D - 0.1) == Selection::Outlier);
  BOOST_CHECK(selectSingle<3>(selector, outlier3D + 0.1) ==
              Selection::Rejected);
}

/// If the chi2 cut depends on the dimension, the candidate with the lowest chi2
/// can fail its cut while a candidate with a higher chi2 passes. The selector
/// must return the candidate that passes.
BOOST_AUTO_TEST_CASE(Chi2CutPerDimensionMixedDimensions) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  auto surface = makeSurface();

  MeasurementSelectorCuts cuts;
  cuts.chi2CutOff = {15};
  cuts.chi2CutOffOutlier = {25};
  cuts.chi2CutPerDimension = true;
  MeasurementSelector selector{cuts};

  // The 1D candidate fails the 1D cut of 15. The 3D candidate has a higher
  // chi2 but passes the 3D cut of 20.96. Check both orders of the candidates.
  for (bool oneDimFirst : {true, false}) {
    BOOST_TEST_CONTEXT("oneDimFirst=" << oneDimFirst) {
      VectorMultiTrajectory traj;
      auto ts1D = addCandidateWithChi2<1>(traj, surface, 16);
      auto ts3D = addCandidateWithChi2<3>(traj, surface, 19);
      std::vector<VectorMultiTrajectory::TrackStateProxy> candidates;
      if (oneDimFirst) {
        candidates = {ts1D, ts3D};
      } else {
        candidates = {ts3D, ts1D};
      }

      bool isOutlier = false;
      auto result = selector.select<VectorMultiTrajectory>(
          surface->boundParametrization(gctx), candidates, isOutlier, *logger);
      BOOST_REQUIRE(result.ok());
      auto [begin, end] = *result;
      BOOST_REQUIRE_EQUAL(end - begin, 1);
      BOOST_CHECK(!isOutlier);
      BOOST_CHECK_EQUAL(begin->calibratedSize(), 3u);
    }
  }
}

/// A NaN or negative chi2 cut has no p-value, so the per dimension conversion
/// rejects it.
BOOST_AUTO_TEST_CASE(Chi2CutPerDimensionInvalidCut) {
  for (double cut : {std::numeric_limits<double>::quiet_NaN(), -1.}) {
    MeasurementSelectorCuts cuts;
    cuts.chi2CutOff = {cut};
    cuts.chi2CutPerDimension = true;
    BOOST_CHECK_THROW(MeasurementSelector{cuts}, std::invalid_argument);

    cuts.chi2CutPerDimension = false;
    BOOST_CHECK_NO_THROW(MeasurementSelector{cuts});
  }
}

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
  auto result = selector.select<VectorMultiTrajectory>(
      surface->boundParametrization(gctx), candidates, isOutlier, *logger);

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
  auto result = selector.select<VectorMultiTrajectory>(
      surface->boundParametrization(gctx), candidates, isOutlier, *logger);

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
        surface->boundParametrization(gctx), candidates, isOutlier, *logger);
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
      surface->boundParametrization(gctx), candidates, isOutlier, *logger);
  BOOST_REQUIRE(result.ok());
  BOOST_REQUIRE_EQUAL(result->second - result->first, 1);
  BOOST_CHECK_EQUAL(result->first->index(), valid.index());
  BOOST_CHECK_SMALL(valid.chi2() - 2.5, 1e-6);
  BOOST_CHECK(!std::isfinite(invalid.chi2()));
  BOOST_CHECK(isOutlier);

  // Position residuals must not wrap like direction coordinates.
  valid.calibrated<2>()[1] = 2. * std::numbers::pi + .1;
  const auto rejected = selector.select<VectorMultiTrajectory>(
      surface->boundParametrization(gctx), candidates, isOutlier, *logger);
  BOOST_REQUIRE(rejected.ok());
  BOOST_CHECK(rejected->first == rejected->second);
  BOOST_CHECK(!isOutlier);
}

/// Candidates sharing one finite prediction are all evaluated, even though
/// the prediction is only checked once.
BOOST_AUTO_TEST_CASE(SharedFinitePrediction) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  auto surface = makeSurface();

  VectorMultiTrajectory traj;
  auto first = addCandidate(traj, surface, false);
  setMeasurement(first, .2);
  auto second = addSharedCandidate(traj, first, .1);
  auto third = addSharedCandidate(traj, first, .3);
  BOOST_REQUIRE_EQUAL(second.predicted().data(), first.predicted().data());
  BOOST_REQUIRE_EQUAL(third.predicted().data(), first.predicted().data());
  std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{first, second,
                                                                 third};

  const auto selector = makeSelector(3);
  bool isOutlier = true;
  const auto result = selector.select<VectorMultiTrajectory>(
      surface->boundParametrization(gctx), candidates, isOutlier, *logger);

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(!isOutlier);
  BOOST_REQUIRE_EQUAL(result->second - result->first, 3);
  // Sorted by chi2.
  BOOST_CHECK_EQUAL(candidates[0].index(), second.index());
  BOOST_CHECK_EQUAL(candidates[1].index(), first.index());
  BOOST_CHECK_EQUAL(candidates[2].index(), third.index());
  BOOST_CHECK_SMALL(second.chi2() - .5, 1e-9);
  BOOST_CHECK_SMALL(first.chi2() - 2., 1e-9);
  BOOST_CHECK_SMALL(third.chi2() - 4.5, 1e-9);
}

/// Candidates sharing one non-finite prediction are all rejected.
BOOST_AUTO_TEST_CASE(SharedNonFinitePrediction) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  auto surface = makeSurface();

  for (double bad : {kNaN, std::numeric_limits<double>::infinity()}) {
    VectorMultiTrajectory traj;
    auto first = addCandidate(traj, surface, false);
    first.predicted()[eBoundDir0] = bad;
    auto second = addSharedCandidate(traj, first, .1);
    auto third = addSharedCandidate(traj, first, .2);
    std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{
        first, second, third};

    const auto selector = makeSelector(3);
    bool isOutlier = true;
    const auto result = selector.select<VectorMultiTrajectory>(
        surface->boundParametrization(gctx), candidates, isOutlier, *logger);

    BOOST_REQUIRE(result.ok());
    BOOST_CHECK(!isOutlier);
    BOOST_CHECK(result->first == result->second);
    for (const auto& candidate : candidates) {
      BOOST_CHECK(std::isnan(candidate.chi2()));
    }
  }
}

/// Candidates with separate predictions are checked individually, also when
/// shared and separate predictions alternate.
BOOST_AUTO_TEST_CASE(SeparatePredictionsOneNonFinite) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  auto surface = makeSurface();

  VectorMultiTrajectory traj;
  auto a = addCandidate(traj, surface, false);
  setMeasurement(a, .1);
  auto b = addCandidate(traj, surface, true);
  auto c = addSharedCandidate(traj, a, .2);
  auto d = addSharedCandidate(traj, b, 0.);
  auto e = addCandidate(traj, surface, false);
  setMeasurement(e, .3);
  std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{a, b, c, d, e};

  const auto selector = makeSelector(5);
  bool isOutlier = true;
  const auto result = selector.select<VectorMultiTrajectory>(
      surface->boundParametrization(gctx), candidates, isOutlier, *logger);

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(!isOutlier);
  BOOST_REQUIRE_EQUAL(result->second - result->first, 3);
  BOOST_CHECK_EQUAL(candidates[0].index(), a.index());
  BOOST_CHECK_EQUAL(candidates[1].index(), c.index());
  BOOST_CHECK_EQUAL(candidates[2].index(), e.index());
  BOOST_CHECK_SMALL(a.chi2() - .5, 1e-9);
  BOOST_CHECK_SMALL(c.chi2() - 2., 1e-9);
  BOOST_CHECK_SMALL(e.chi2() - 4.5, 1e-9);
  BOOST_CHECK(std::isnan(b.chi2()));
  BOOST_CHECK(std::isnan(d.chi2()));

  // A non-finite prediction ahead of the first finite one.
  std::vector<VectorMultiTrajectory::TrackStateProxy> reordered{b, d, e, a};
  const auto second = selector.select<VectorMultiTrajectory>(
      surface->boundParametrization(gctx), reordered, isOutlier, *logger);
  BOOST_REQUIRE(second.ok());
  BOOST_REQUIRE_EQUAL(second->second - second->first, 2);
  BOOST_CHECK_EQUAL(reordered[0].index(), a.index());
  BOOST_CHECK_EQUAL(reordered[1].index(), e.index());
  BOOST_CHECK(std::isnan(b.chi2()));
  BOOST_CHECK(std::isnan(d.chi2()));
}

/// A non-finite measurement is rejected per candidate, also when its
/// prediction is shared with candidates that are selected.
BOOST_AUTO_TEST_CASE(NonFiniteCalibrated) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  auto surface = makeSurface();

  VectorMultiTrajectory traj;
  auto first = addCandidate(traj, surface, false);
  setMeasurement(first, kNaN);
  auto second = addSharedCandidate(traj, first, .1);
  auto third =
      addSharedCandidate(traj, first, std::numeric_limits<double>::infinity());
  auto fourth = addSharedCandidate(traj, first, .2);
  std::vector<VectorMultiTrajectory::TrackStateProxy> candidates{first, second,
                                                                 third, fourth};

  const auto selector = makeSelector(4);
  bool isOutlier = true;
  const auto result = selector.select<VectorMultiTrajectory>(
      surface->boundParametrization(gctx), candidates, isOutlier, *logger);

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK(!isOutlier);
  BOOST_REQUIRE_EQUAL(result->second - result->first, 2);
  BOOST_CHECK_EQUAL(candidates[0].index(), second.index());
  BOOST_CHECK_EQUAL(candidates[1].index(), fourth.index());
  BOOST_CHECK_SMALL(second.chi2() - .5, 1e-9);
  BOOST_CHECK_SMALL(fourth.chi2() - 2., 1e-9);
  BOOST_CHECK(std::isnan(first.chi2()));
  BOOST_CHECK(std::isnan(third.chi2()));
}

namespace {

struct TestSourceLink {
  double loc0 = 0.;
};

void calibrateTestSourceLink(const GeometryContext& /*gctx*/,
                             const CalibrationContext& /*cctx*/,
                             const SourceLink& sourceLink,
                             VectorMultiTrajectory::TrackStateProxy ts) {
  setMeasurement(ts, sourceLink.get<TestSourceLink>().loc0);
}

/// Forwards to the measurement selector and records the parametrization the
/// track state creator passes.
struct RecordingSelector {
  MeasurementSelector selector = makeSelector(3);
  mutable std::optional<BoundParametrization> parametrization;
  mutable std::size_t calls = 0;

  Result<
      std::pair<std::vector<VectorMultiTrajectory::TrackStateProxy>::iterator,
                std::vector<VectorMultiTrajectory::TrackStateProxy>::iterator>>
  select(const BoundParametrization& passed,
         std::vector<VectorMultiTrajectory::TrackStateProxy>& candidates,
         bool& isOutlier, const Logger& logger) const {
    parametrization = passed;
    ++calls;
    return selector.select<VectorMultiTrajectory>(passed, candidates, isOutlier,
                                                  logger);
  }
};

}  // namespace

/// The track state creator resolves the parametrization of the candidates'
/// surface once and passes it to the selector.
BOOST_AUTO_TEST_CASE(TrackStateCreatorPassesSurfaceParametrization) {
  auto logger = getDefaultLogger("MeasurementSelector", Logging::INFO);
  const auto surface = makeSurface();
  const CalibrationContext cctx;

  using TrackContainer =
      Acts::TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                           std::shared_ptr>;
  using Creator = TrackStateCreator<std::vector<SourceLink>::const_iterator,
                                    TrackContainer>;
  Creator creator;
  creator.calibrator.connect<&calibrateTestSourceLink>();
  RecordingSelector recorder;
  creator.measurementSelector.connect<&RecordingSelector::select>(&recorder);

  BoundVector parameters = BoundVector::Zero();
  parameters[eBoundTheta] = std::numbers::pi / 2.;
  parameters[eBoundQOverP] = 1.;
  const Creator::BoundState boundState{
      BoundTrackParameters(surface, parameters, BoundMatrix::Identity() * 1e-2,
                           ParticleHypothesis::pion()),
      BoundMatrix::Identity(), 1.};
  const std::vector<SourceLink> sourceLinks{SourceLink{TestSourceLink{.2}},
                                            SourceLink{TestSourceLink{.1}},
                                            SourceLink{TestSourceLink{10.}}};

  VectorMultiTrajectory trajectory;
  std::vector<Creator::TrackStateProxy> candidates;
  const auto result = creator.createSourceLinkTrackStates(
      gctx, cctx, *surface, boundState, sourceLinks.begin(), sourceLinks.end(),
      kTrackIndexInvalid, candidates, trajectory, *logger);

  BOOST_REQUIRE(result.ok());
  BOOST_CHECK_EQUAL(recorder.calls, 1u);
  BOOST_REQUIRE(recorder.parametrization.has_value());
  BOOST_CHECK(*recorder.parametrization == surface->boundParametrization(gctx));
  // chi2 .5 and 2 pass the cut of 15, the third candidate does not.
  BOOST_REQUIRE_EQUAL(result->size(), 2u);
  BOOST_CHECK_SMALL(trajectory.getTrackState(result->at(0)).chi2() - .5, 1e-9);
  BOOST_CHECK_SMALL(trajectory.getTrackState(result->at(1)).chi2() - 2., 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
