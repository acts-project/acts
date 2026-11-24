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
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(MeasurementSelectorTests)

namespace {

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
  ts.predicted()[eBoundTheta] = std::numbers::pi / 2.;
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
  auto result =
      selector.select<VectorMultiTrajectory>(candidates, isOutlier, *logger);
  BOOST_REQUIRE(result.ok());
  auto [begin, end] = *result;
  BOOST_REQUIRE_LE(end - begin, 1);
  if (begin == end) {
    BOOST_CHECK(!isOutlier);
    return Selection::Rejected;
  }
  return isOutlier ? Selection::Outlier : Selection::Measurement;
}

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
      auto result = selector.select<VectorMultiTrajectory>(candidates,
                                                           isOutlier, *logger);
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
  auto result =
      selector.select<VectorMultiTrajectory>(candidates, isOutlier, *logger);

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
  auto result =
      selector.select<VectorMultiTrajectory>(candidates, isOutlier, *logger);

  BOOST_REQUIRE(result.ok());
  auto [begin, end] = *result;
  BOOST_CHECK_GE(begin - candidates.begin(), 0);
  BOOST_REQUIRE_EQUAL(end - begin, 1);
  BOOST_CHECK(std::isfinite(begin->chi2()));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
