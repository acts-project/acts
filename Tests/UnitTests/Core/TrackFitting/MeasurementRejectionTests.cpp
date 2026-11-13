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
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"
#include "Acts/TrackFitting/MeasurementRejection.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {

using namespace Acts;

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

using Container = TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                 detail::ValueHolder>;

constexpr std::size_t kNoOutlier = static_cast<std::size_t>(-1);

enum class Kind { Measurement2, Measurement1, Hole };

/// One node of a linear-Gaussian chain: transport from the previous node, an
/// affine offset as from the mean energy loss, process noise, and a
/// measurement of (loc0, loc1) or loc0.
struct Node {
  BoundMatrix F = BoundMatrix::Identity();
  BoundVector d = BoundVector::Zero();
  BoundMatrix Q = BoundMatrix::Zero();
  Kind kind = Kind::Measurement2;
  Vector2 y = Vector2::Zero();
  SquareMatrix2 V = SquareMatrix2::Identity();
};

std::vector<Node> makeChain() {
  const std::array<Kind, 9> kinds = {
      Kind::Measurement2, Kind::Hole,         Kind::Measurement2,
      Kind::Measurement1, Kind::Measurement2, Kind::Hole,
      Kind::Measurement2, Kind::Measurement1, Kind::Measurement2};
  std::vector<Node> nodes;
  for (std::size_t i = 0; i < kinds.size(); ++i) {
    const double s = 0.7 + 0.1 * static_cast<double>(i % 4);
    Node n;
    n.F(eBoundLoc0, eBoundPhi) = 0.3 * s;
    n.F(eBoundLoc1, eBoundTheta) = 0.2 * s;
    n.F(eBoundPhi, eBoundQOverP) = 0.1 * s;
    n.F(eBoundLoc0, eBoundQOverP) = 0.05 * s;
    n.F(eBoundTime, eBoundLoc1) = 0.1 * s;
    n.d[eBoundQOverP] = -1e-3 * s;
    BoundVector q;
    q << 1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-3;
    n.Q = (q * s).asDiagonal();
    n.kind = kinds[i];
    n.y = Vector2(0.1 * static_cast<double>(i) - 0.2,
                  0.05 * static_cast<double>(i * i % 5) - 0.1);
    n.V = Vector2(4e-3 * s, 9e-3 / s).asDiagonal();
    nodes.push_back(n);
  }
  return nodes;
}

const BoundVector kPriorMean = [] {
  BoundVector m;
  m << 0.1, -0.2, 1.3, 0.7, 0.05, 0.0;
  return m;
}();

const BoundMatrix kPriorCov = [] {
  BoundVector d;
  d << 0.04, 0.09, 0.01, 0.01, 1e-4, 0.02;
  return BoundMatrix(d.asDiagonal());
}();

/// Filter the chain the way the CKF stores it, with @p outlier as outlier,
/// and smooth it
template <typename smoother_t = MbfSmoother>
Container::TrackProxy buildTrack(Container& tc, const std::vector<Node>& nodes,
                                 std::size_t outlier = kNoOutlier) {
  auto track = tc.makeTrack();
  BoundVector mean = kPriorMean;
  BoundMatrix cov = kPriorCov;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const Node& node = nodes[i];
    auto ts = track.appendTrackState(TrackStatePropMask::All);
    ts.jacobian() = i == 0 ? BoundMatrix::Identity() : node.F;
    if (i > 0) {
      mean = node.F * mean + node.d;
      cov = node.F * cov * node.F.transpose() + node.Q;
    }
    ts.predicted() = mean;
    ts.predictedCovariance() = cov;
    if (node.kind == Kind::Hole) {
      ts.typeFlags().setIsHole();
      ts.shareFrom(TrackStatePropMask::Predicted, TrackStatePropMask::Filtered);
      continue;
    }
    const std::size_t dim = node.kind == Kind::Measurement2 ? 2 : 1;
    ts.allocateCalibrated(dim);
    if (dim == 2) {
      ts.calibrated<2>() = node.y;
      ts.calibratedCovariance<2>() = node.V;
      ts.setProjectorSubspaceIndices(std::array{eBoundLoc0, eBoundLoc1});
    } else {
      ts.calibrated<1>() = node.y.head<1>();
      ts.calibratedCovariance<1>() = node.V.topLeftCorner<1, 1>();
      ts.setProjectorSubspaceIndices(std::array{eBoundLoc0});
    }
    const DynamicMatrix h = BoundMatrix::Identity().topRows(dim);
    const DynamicMatrix s =
        h * cov * h.transpose() + node.V.topLeftCorner(dim, dim);
    const DynamicVector r = node.y.head(dim) - h * mean;
    ts.chi2() = static_cast<float>(r.dot(s.inverse() * r));
    if (i == outlier) {
      ts.typeFlags().setIsOutlier();
      ts.shareFrom(TrackStatePropMask::Predicted, TrackStatePropMask::Filtered);
      track.nOutliers()++;
      continue;
    }
    ts.typeFlags().setIsMeasurement();
    const DynamicMatrix k = cov * h.transpose() * s.inverse();
    mean = mean + k * r;
    cov = cov - k * h * cov;
    ts.filtered() = mean;
    ts.filteredCovariance() = cov;
    track.nMeasurements()++;
    track.nDoF() += dim;
    track.chi2() += ts.chi2();
  }
  BOOST_REQUIRE(
      smoother_t{}(gctx, tc.trackStateContainer(), track.tipIndex()).ok());
  return track;
}

template <typename smoother_t>
void checkRejection(std::size_t k) {
  const auto nodes = makeChain();
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto track = buildTrack<smoother_t>(tc, nodes);
  auto expected = buildTrack(tc, nodes, k);

  auto state = track.trackStatesReversed().begin();
  for (std::size_t i = nodes.size() - 1; i > k; --i) {
    ++state;
  }
  BOOST_REQUIRE(rejectMeasurement(track, *state).ok());
  BOOST_CHECK((*state).typeFlags().isOutlier());
  CHECK_CLOSE_ABS((*state).filtered(), (*state).predicted(), 0.);

  BOOST_CHECK_EQUAL(track.nMeasurements(), expected.nMeasurements());
  BOOST_CHECK_EQUAL(track.nOutliers(), expected.nOutliers());
  BOOST_CHECK_EQUAL(track.nDoF(), expected.nDoF());
  CHECK_CLOSE_REL(track.chi2(), expected.chi2(), 1e-5);

  auto expectedStates = expected.trackStatesReversed();
  auto e = expectedStates.begin();
  for (const auto& ts : track.trackStatesReversed()) {
    CHECK_CLOSE_ABS(ts.smoothed(), (*e).smoothed(), 1e-10);
    CHECK_CLOSE_ABS(ts.smoothedCovariance(), (*e).smoothedCovariance(), 1e-10);
    ++e;
  }
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

using Smoothers = std::tuple<MbfSmoother, GainMatrixSmoother>;

BOOST_AUTO_TEST_CASE_TEMPLATE(RejectionMatchesFitWithOutlier, smoother_t,
                              Smoothers) {
  // first, last, one- and two-dimensional, next to holes
  for (std::size_t k : {0u, 2u, 3u, 4u, 6u, 7u, 8u}) {
    BOOST_TEST_CONTEXT("rejected node " << k) {
      checkRejection<smoother_t>(k);
    }
  }
}

BOOST_AUTO_TEST_CASE(RejectionNeedsMeasurementOnTrack) {
  const auto nodes = makeChain();
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto track = buildTrack(tc, nodes);
  auto other = buildTrack(tc, nodes);
  // the tip is a measurement, the state before it a hole
  auto hole = *std::next(track.trackStatesReversed().begin(), 3);
  BOOST_REQUIRE(hole.typeFlags().isHole());
  BOOST_CHECK_THROW((void)rejectMeasurement(track, hole),
                    std::invalid_argument);
  BOOST_CHECK_THROW((void)rejectMeasurement(track, other.outermostTrackState()),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RejectionRejectsNonFiniteInputs) {
  for (const double value : {std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity()}) {
    for (const bool covariance : {false, true}) {
      Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
      auto nodes = makeChain();
      nodes.resize(1);
      auto track = buildTrack(tc, nodes);
      auto state = track.outermostTrackState();
      if (covariance) {
        state.calibratedCovariance<2>()(0, 0) = value;
      } else {
        state.calibrated<2>()[0] = value;
      }
      BOOST_CHECK(!rejectMeasurement(track, state).ok());
      BOOST_CHECK(state.smoothed().allFinite());
      BOOST_CHECK(state.typeFlags().isMeasurement());
    }
  }
}

BOOST_AUTO_TEST_CASE(RejectionRejectsStateFromAnotherTrajectory) {
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  Container other{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto track = buildTrack(tc, makeChain());
  auto foreign = buildTrack(other, makeChain()).outermostTrackState();
  BOOST_CHECK_THROW((void)rejectMeasurement(track, foreign),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(RejectionWithPermutedProjectors, smoother_t,
                              Smoothers) {
  // Exercise all update ranks and non-contiguous, reordered projectors.
  for (std::size_t offset = 0; offset < eBoundSize; ++offset) {
    Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
    auto build = [&](bool omit) {
      auto track = tc.makeTrack();
      BoundVector mean = kPriorMean;
      BoundMatrix cov = kPriorCov + 0.001 * BoundMatrix::Ones();
      for (std::size_t i = 0; i < 7; ++i) {
        auto ts = track.appendTrackState(TrackStatePropMask::All);
        BoundMatrix f = BoundMatrix::Identity();
        f(eBoundLoc0, eBoundPhi) = 0.03;
        f(eBoundLoc1, eBoundTheta) = 0.02;
        mean = f * mean;
        cov = f * cov * f.transpose() + 0.001 * BoundMatrix::Identity();
        ts.jacobian() = f;
        ts.predicted() = mean;
        ts.predictedCovariance() = cov;
        const std::size_t dim = 1 + (i + offset) % eBoundSize;
        std::vector<BoundIndices> indices;
        DynamicMatrix h = DynamicMatrix::Zero(dim, eBoundSize);
        for (std::size_t j = 0; j < dim; ++j) {
          indices.push_back(
              static_cast<BoundIndices>((5 * j + i) % eBoundSize));
          h(j, indices.back()) = 1;
        }
        ts.allocateCalibrated(dim);
        ts.setProjectorSubspaceIndices(indices);
        const DynamicMatrix v = 0.01 * DynamicMatrix::Identity(dim, dim);
        const DynamicVector y =
            h * kPriorMean + DynamicVector::Constant(dim, 0.01 * i);
        ts.effectiveCalibrated() = y;
        ts.effectiveCalibratedCovariance() = v;
        if (omit && i == 3) {
          ts.typeFlags().setIsOutlier();
          ts.shareFrom(TrackStatePropMask::Predicted,
                       TrackStatePropMask::Filtered);
          track.nOutliers()++;
          continue;
        }
        ts.typeFlags().setIsMeasurement();
        const DynamicMatrix sInv = (h * cov * h.transpose() + v).inverse();
        const DynamicMatrix k = cov * h.transpose() * sInv;
        const DynamicVector residual = y - h * mean;
        mean += k * residual;
        cov -= k * h * cov;
        ts.filtered() = mean;
        ts.filteredCovariance() = cov;
        ts.chi2() = static_cast<float>(residual.dot(sInv * residual));
        track.chi2() += ts.chi2();
        track.nDoF() += dim;
        track.nMeasurements()++;
      }
      BOOST_REQUIRE(
          smoother_t{}(gctx, tc.trackStateContainer(), track.tipIndex()).ok());
      return track;
    };
    auto track = build(false);
    auto expected = build(true);
    auto state = *std::next(track.trackStatesReversed().begin(), 3);
    BOOST_REQUIRE(rejectMeasurement(track, state).ok());
    CHECK_CLOSE_ABS(track.chi2(), expected.chi2(), 1e-5);
    auto e = expected.trackStatesReversed().begin();
    for (const auto& ts : track.trackStatesReversed()) {
      CHECK_CLOSE_ABS(ts.smoothed(), (*e).smoothed(), 1e-10);
      CHECK_CLOSE_ABS(ts.smoothedCovariance(), (*e).smoothedCovariance(),
                      1e-10);
      ++e;
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
