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
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"
#include "Acts/TrackFitting/SmootherAdjoint.hpp"
#include "Acts/TrackFitting/TrackletJoiner.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <vector>

#include <Eigen/Dense>

namespace {

using namespace Acts;

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

constexpr std::size_t kDim = eBoundSize;
const std::array<BoundIndices, 2> kProjector{eBoundLoc0, eBoundLoc1};

/// One node of a linear-Gaussian model: transport from its parent, process
/// noise, and an optional two-dimensional measurement of (loc0, loc1).
struct Node {
  int parent = -1;  // -1: the node carries the prior
  BoundMatrix F = BoundMatrix::Identity();
  BoundMatrix Q = BoundMatrix::Zero();
  bool hasMeasurement = false;
  Vector2 y = Vector2::Zero();
  SquareMatrix2 V = SquareMatrix2::Identity();
};

BoundMatrix makeJacobian(double scale) {
  BoundMatrix f = BoundMatrix::Identity();
  f(eBoundLoc0, eBoundPhi) = 0.3 * scale;
  f(eBoundLoc1, eBoundTheta) = 0.2 * scale;
  f(eBoundPhi, eBoundQOverP) = 0.1 * scale;
  f(eBoundLoc0, eBoundQOverP) = 0.05 * scale;
  f(eBoundTime, eBoundLoc1) = 0.1 * scale;
  return f;
}

BoundMatrix makeNoise(double scale) {
  BoundVector d;
  d << 1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-3;
  return (d * scale).asDiagonal();
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

Eigen::Matrix<double, 2, kDim> projector() {
  Eigen::Matrix<double, 2, kDim> h = Eigen::Matrix<double, 2, kDim>::Zero();
  h(0, eBoundLoc0) = 1;
  h(1, eBoundLoc1) = 1;
  return h;
}

/// Variants of the model
enum class Variant {
  /// nodes 0 to 3 are the first tracklet in filter order, node 0 carries the
  /// prior and is the shared state; nodes 4 to 6 are the second tracklet,
  /// filtered from node 0. Nodes 1 and 5 are holes.
  Base,
  /// Base plus node 7, a hole after the last measurement of the first
  /// tracklet. It has no smoothed component, as after `smoothTrack`.
  TrailingHole,
  /// Base plus node 7, a hole before the first measurement. It carries the
  /// prior, node 0 is transported from it, and node 0 stays the shared state.
  StemHole,
};

/// The model
std::vector<Node> makeModel(Variant variant) {
  std::vector<Node> nodes(7);
  const SquareMatrix2 v = Vector2(4e-3, 9e-3).asDiagonal();

  nodes[0] = {-1,   BoundMatrix::Identity(), BoundMatrix::Zero(),
              true, Vector2(0.12, -0.23),    v};
  nodes[1] = {0, makeJacobian(1.0), makeNoise(1.0), false, Vector2::Zero(), v};
  nodes[2] = {1,    makeJacobian(0.8),   makeNoise(1.5),
              true, Vector2(0.55, -0.1), v};
  nodes[3] = {2,    makeJacobian(1.2),  makeNoise(0.7),
              true, Vector2(0.9, 0.05), v};
  nodes[4] = {0,    makeJacobian(-1.1),    makeNoise(1.2),
              true, Vector2(-0.25, -0.35), v};
  nodes[5] = {4, makeJacobian(-0.9), makeNoise(0.9), false, Vector2::Zero(), v};
  nodes[6] = {5,    makeJacobian(-1.3),  makeNoise(1.1),
              true, Vector2(-0.7, -0.5), v};

  if (variant == Variant::TrailingHole) {
    nodes.push_back(
        {3, makeJacobian(0.9), makeNoise(1.3), false, Vector2::Zero(), v});
  } else if (variant == Variant::StemHole) {
    // the stem hole carries the prior, node 0 is transported from it
    nodes.push_back({-1, BoundMatrix::Identity(), BoundMatrix::Zero(), false,
                     Vector2::Zero(), v});
    nodes[0].parent = 7;
    nodes[0].F = makeJacobian(0.7);
    nodes[0].Q = makeNoise(0.8);
  }
  return nodes;
}

/// Dense joint posterior of all nodes given all measurements and the prior.
std::pair<Eigen::VectorXd, Eigen::MatrixXd> densePosterior(
    const std::vector<Node>& nodes) {
  const std::size_t n = nodes.size();
  const std::size_t dim = n * kDim;

  Eigen::MatrixXd a = Eigen::MatrixXd::Zero(dim, dim);
  Eigen::MatrixXd w = Eigen::MatrixXd::Zero(dim, dim);
  Eigen::VectorXd c = Eigen::VectorXd::Zero(dim);
  for (std::size_t i = 0; i < n; ++i) {
    const Node& node = nodes[i];
    if (node.parent < 0) {
      w.block<kDim, kDim>(i * kDim, i * kDim) = kPriorCov;
      c.segment<kDim>(i * kDim) = kPriorMean;
    } else {
      a.block<kDim, kDim>(i * kDim, node.parent * kDim) = node.F;
      w.block<kDim, kDim>(i * kDim, i * kDim) = node.Q;
    }
  }

  const Eigen::MatrixXd g = (Eigen::MatrixXd::Identity(dim, dim) - a).inverse();
  const Eigen::VectorXd priorMean = g * c;
  const Eigen::MatrixXd priorCov = g * w * g.transpose();

  Eigen::MatrixXd information = priorCov.inverse();
  Eigen::VectorXd informationVector = information * priorMean;
  const Eigen::Matrix<double, 2, kDim> h = projector();
  for (std::size_t i = 0; i < n; ++i) {
    const Node& node = nodes[i];
    if (!node.hasMeasurement) {
      continue;
    }
    const SquareMatrix2 vInv = node.V.inverse();
    information.block<kDim, kDim>(i * kDim, i * kDim) +=
        h.transpose() * vInv * h;
    informationVector.segment<kDim>(i * kDim) += h.transpose() * vInv * node.y;
  }

  const Eigen::MatrixXd cov = information.inverse();
  const Eigen::VectorXd mean = cov * informationVector;
  return {mean, cov};
}

using Container = TrackContainer<VectorTrackContainer, VectorMultiTrajectory,
                                 detail::ValueHolder>;

/// Append one filtered state to a track, the way the CKF stores it.
std::size_t appendState(Container::TrackProxy& track, const Node& node,
                        const BoundVector& prevMean, const BoundMatrix& prevCov,
                        TrackStatePropMask mask = TrackStatePropMask::All) {
  auto ts = track.appendTrackState(mask);
  ts.jacobian() = node.F;

  const BoundVector predicted = node.parent < 0 ? prevMean : node.F * prevMean;
  const BoundMatrix predictedCov =
      node.parent < 0 ? prevCov
                      : node.F * prevCov * node.F.transpose() + node.Q;
  ts.predicted() = predicted;
  ts.predictedCovariance() = predictedCov;

  if (node.hasMeasurement) {
    ts.typeFlags().setIsMeasurement();
    ts.allocateCalibrated(2);
    ts.calibrated<2>() = node.y;
    ts.calibratedCovariance<2>() = node.V;
    ts.setProjectorSubspaceIndices(kProjector);

    const Eigen::Matrix<double, 2, kDim> h = projector();
    const SquareMatrix2 s = h * predictedCov * h.transpose() + node.V;
    const Eigen::Matrix<double, kDim, 2> k =
        predictedCov * h.transpose() * s.inverse();
    ts.filtered() = predicted + k * (node.y - h * predicted);
    ts.filteredCovariance() = (BoundMatrix::Identity() - k * h) * predictedCov;
  } else {
    ts.typeFlags().setIsHole();
    ts.filtered() = predicted;
    ts.filteredCovariance() = predictedCov;
  }
  return ts.index();
}

struct Tracklets {
  Container::TrackProxy first;
  Container::TrackStateProxy shared;
  Container::TrackProxy second;
  /// node of each track state index
  std::map<std::size_t, std::size_t> nodeOfState;
  /// states that must stay without smoothed parameters
  std::vector<std::size_t> unsmoothed;
};

/// Build the first tracklet, smooth it from its last measurement as
/// `smoothTrack` does, and build the second tracklet from the smoothed
/// shared state, which is the first measurement.
Tracklets buildTracklets(Container& tc, const std::vector<Node>& nodes,
                         Variant variant) {
  std::map<std::size_t, std::size_t> nodeOfState;
  std::vector<std::size_t> unsmoothed;

  auto first = tc.makeTrack();
  BoundVector mean = kPriorMean;
  BoundMatrix cov = kPriorCov;
  std::vector<std::size_t> firstNodes = {0, 1, 2, 3};
  if (variant == Variant::TrailingHole) {
    firstNodes.push_back(7);
  } else if (variant == Variant::StemHole) {
    firstNodes.insert(firstNodes.begin(), 7);
  }
  std::size_t lastMeasurement = kTrackIndexInvalid;
  std::size_t sharedIndex = kTrackIndexInvalid;
  for (std::size_t node : firstNodes) {
    // a trailing hole is created without a smoothed component, as the CKF
    // does for every state
    const bool trailing = (variant == Variant::TrailingHole && node == 7);
    const std::size_t index = appendState(
        first, nodes[node], mean, cov,
        trailing ? (TrackStatePropMask::Predicted |
                    TrackStatePropMask::Filtered | TrackStatePropMask::Jacobian)
                 : TrackStatePropMask::All);
    nodeOfState[index] = node;
    auto ts = first.outermostTrackState();
    mean = ts.filtered();
    cov = ts.filteredCovariance();
    if (nodes[node].hasMeasurement) {
      lastMeasurement = index;
      if (sharedIndex == kTrackIndexInvalid) {
        sharedIndex = index;
      }
    }
    if (trailing) {
      unsmoothed.push_back(index);
    }
  }
  first.linkForward();
  BOOST_REQUIRE(
      MbfSmoother()(gctx, tc.trackStateContainer(), lastMeasurement).ok());

  auto shared = tc.trackStateContainer().getTrackState(sharedIndex);
  BOOST_REQUIRE(shared.hasSmoothed());
  BOOST_REQUIRE(shared.typeFlags().isMeasurement());

  auto second = tc.makeTrack();
  mean = shared.smoothed();
  cov = shared.smoothedCovariance();
  for (std::size_t node = 4; node < 7; ++node) {
    const std::size_t index = appendState(second, nodes[node], mean, cov);
    nodeOfState[index] = node;
    auto ts = second.outermostTrackState();
    mean = ts.filtered();
    cov = ts.filteredCovariance();
  }

  return {first, shared, second, nodeOfState, unsmoothed};
}

/// Check every smoothed state of a joined track against the dense joint
/// posterior. States listed as unsmoothed must have no smoothed component.
template <typename track_state_proxy_t>
void checkStateAgainstDense(const track_state_proxy_t& ts, std::size_t node,
                            const Eigen::VectorXd& mean,
                            const Eigen::MatrixXd& cov, double tol) {
  BOOST_REQUIRE(ts.hasSmoothed());
  const BoundVector expMean = mean.segment<kDim>(node * kDim);
  const BoundMatrix expCov = cov.block<kDim, kDim>(node * kDim, node * kDim);
  CHECK_CLOSE_ABS(ts.smoothed(), expMean, tol);
  CHECK_CLOSE_ABS(ts.smoothedCovariance(), expCov, tol);
  // the covariance is symmetric
  const BoundMatrix asym =
      ts.smoothedCovariance() - ts.smoothedCovariance().transpose();
  BOOST_CHECK_LT(asym.cwiseAbs().maxCoeff(), 1e-15);
}

void checkJoinedAgainstDense(const Container::TrackProxy& joined,
                             const std::map<std::size_t, std::size_t>& nodeOf,
                             const std::vector<std::size_t>& unsmoothed,
                             const std::vector<Node>& nodes, double tol,
                             std::size_t expectedStates) {
  const auto [mean, cov] = densePosterior(nodes);
  std::size_t count = 0;
  for (const auto& ts : joined.trackStatesReversed()) {
    ++count;
    if (std::find(unsmoothed.begin(), unsmoothed.end(), ts.index()) !=
        unsmoothed.end()) {
      BOOST_CHECK(!ts.hasSmoothed());
      continue;
    }
    BOOST_REQUIRE(nodeOf.contains(ts.index()));
    checkStateAgainstDense(ts, nodeOf.at(ts.index()), mean, cov, tol);
  }
  BOOST_CHECK_EQUAL(count, expectedStates);
}

/// Join in place and check against the dense posterior
template <typename smoother_t>
Tracklets joinAndCheck(Container& tc, Variant variant) {
  const auto nodes = makeModel(variant);
  auto tracklets = buildTracklets(tc, nodes, variant);
  TrackletJoiner joiner;
  BOOST_REQUIRE(joiner
                    .joinInPlace(gctx, tracklets.first, tracklets.shared,
                                 tracklets.second, getDummyLogger(),
                                 smoother_t{})
                    .ok());
  // states before the shared state leave the chain in the link
  const std::size_t onChain = variant == Variant::StemHole ? 7 : nodes.size();
  checkJoinedAgainstDense(tracklets.first, tracklets.nodeOfState,
                          tracklets.unsmoothed, nodes, 1e-9, onChain);
  return tracklets;
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(TrackFittingSuite)

using Smoothers = std::tuple<MbfSmoother, GainMatrixSmoother>;

BOOST_AUTO_TEST_CASE_TEMPLATE(CarriedAdjointIsPriorNodeSmoothing, smoother_t,
                              Smoothers) {
  // smoothing with a carried adjoint applied to the node before the segment
  // equals smoothing a chain that holds that node as a state
  const auto nodes = makeModel(Variant::Base);
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto full = tc.makeTrack();
  BoundVector mean = kPriorMean;
  BoundMatrix cov = kPriorCov;
  for (std::size_t i = 0; i < 4; ++i) {
    appendState(full, nodes[i], mean, cov);
    mean = full.outermostTrackState().filtered();
    cov = full.outermostTrackState().filteredCovariance();
  }
  // a copy, smoothed from the tip down to node 2 with the carry
  auto segment = tc.makeTrack();
  segment.copyFrom(full);
  std::vector<std::size_t> segmentIndices;
  for (const auto& ts : segment.trackStatesReversed()) {
    segmentIndices.push_back(ts.index());
  }
  // chain positions from the tip: node 3, 2, 1, 0
  const std::size_t stopIndex = segmentIndices[1];
  const std::size_t priorIndex = segmentIndices[2];

  BOOST_REQUIRE(
      MbfSmoother()(gctx, tc.trackStateContainer(), full.tipIndex()).ok());

  SmootherAdjoint adjoint;
  BOOST_REQUIRE(smoother_t{}
                    .smooth(gctx, tc.trackStateContainer(), segment.tipIndex(),
                            stopIndex, &adjoint)
                    .ok());

  auto priorState = tc.trackStateContainer().getTrackState(priorIndex);
  const auto [priorMean, priorCov] = applySmootherAdjoint(
      priorState.filtered(), priorState.filteredCovariance(), adjoint);

  // the full chain holds the prior node as state at the same position
  std::vector<std::size_t> fullIndices;
  for (const auto& ts : full.trackStatesReversed()) {
    fullIndices.push_back(ts.index());
  }
  auto fullPrior = tc.trackStateContainer().getTrackState(fullIndices[2]);
  CHECK_CLOSE_ABS(priorMean, fullPrior.smoothed(), 1e-10);
  CHECK_CLOSE_ABS(priorCov, fullPrior.smoothedCovariance(), 1e-10);

  // the states down to the stop are smoothed as in the full run
  for (std::size_t i = 0; i < 2; ++i) {
    auto a = tc.trackStateContainer().getTrackState(segmentIndices[i]);
    auto b = tc.trackStateContainer().getTrackState(fullIndices[i]);
    CHECK_CLOSE_ABS(a.smoothed(), b.smoothed(), 1e-10);
    CHECK_CLOSE_ABS(a.smoothedCovariance(), b.smoothedCovariance(), 1e-10);
  }

  // an unreachable stop index is an error
  BOOST_CHECK(!smoother_t{}
                   .smooth(gctx, tc.trackStateContainer(), segment.tipIndex(),
                           full.tipIndex(), nullptr)
                   .ok());
}

BOOST_AUTO_TEST_CASE(CarriedAdjointAgreesBetweenSmoothers) {
  const auto nodes = makeModel(Variant::Base);
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto tracklets = buildTracklets(tc, nodes, Variant::Base);
  auto secondCopy = tc.makeTrack();
  secondCopy.copyFrom(tracklets.second);

  SmootherAdjoint mbf;
  SmootherAdjoint rts;
  BOOST_REQUIRE(MbfSmoother{}
                    .smooth(gctx, tc.trackStateContainer(),
                            tracklets.second.tipIndex(), kTrackIndexInvalid,
                            &mbf)
                    .ok());
  BOOST_REQUIRE(GainMatrixSmoother{}
                    .smooth(gctx, tc.trackStateContainer(),
                            secondCopy.tipIndex(), kTrackIndexInvalid, &rts)
                    .ok());

  CHECK_CLOSE_ABS(mbf.bigLambda, rts.bigLambda, 1e-9);
  CHECK_CLOSE_ABS(mbf.smallLambda, rts.smallLambda, 1e-9);

  const auto [meanA, covA] = applySmootherAdjoint(
      tracklets.shared.smoothed(), tracklets.shared.smoothedCovariance(), mbf);
  const auto [meanB, covB] = applySmootherAdjoint(
      tracklets.shared.smoothed(), tracklets.shared.smoothedCovariance(), rts);
  CHECK_CLOSE_ABS(meanA, meanB, 1e-10);
  CHECK_CLOSE_ABS(covA, covB, 1e-10);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(JoinInPlaceMatchesJointPosterior, smoother_t,
                              Smoothers) {
  const auto nodes = makeModel(Variant::Base);
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto tracklets = buildTracklets(tc, nodes, Variant::Base);
  const BoundMatrix secondStart = nodes[4].F;
  const std::size_t secondInnermost = tracklets.second.tipIndex();

  TrackletJoiner joiner;
  BOOST_REQUIRE(joiner
                    .joinInPlace(gctx, tracklets.first, tracklets.shared,
                                 tracklets.second, getDummyLogger(),
                                 smoother_t{})
                    .ok());

  checkJoinedAgainstDense(tracklets.first, tracklets.nodeOfState,
                          tracklets.unsmoothed, nodes, 1e-9, nodes.size());

  // the link: the chain runs through the shared state into the second
  // tracklet, the stitch jacobian is the inverse transport, the stem is the
  // innermost state of the second tracklet, and the forward links follow
  BOOST_CHECK_EQUAL(tracklets.shared.previous(),
                    tracklets.second.outermostTrackState().index());
  CHECK_CLOSE_ABS(tracklets.shared.jacobian(), secondStart.inverse(), 1e-12);
  BOOST_CHECK_EQUAL(tracklets.first.stemIndex(), secondInnermost);
  std::vector<std::size_t> backward;
  for (const auto& ts : tracklets.first.trackStatesReversed()) {
    backward.push_back(ts.index());
  }
  std::vector<std::size_t> forward;
  for (const auto& ts : tracklets.first.trackStates()) {
    forward.push_back(ts.index());
  }
  BOOST_REQUIRE_EQUAL(backward.size(), 7u);
  BOOST_REQUIRE_EQUAL(forward.size(), 7u);
  for (std::size_t i = 0; i < 7; ++i) {
    BOOST_CHECK_EQUAL(forward[i], backward[6 - i]);
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(JoinSkipsTrailingHoles, smoother_t, Smoothers) {
  // a hole after the last measurement has no smoothed component and stays
  // as it is; every other state matches the joint posterior
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto tracklets = joinAndCheck<smoother_t>(tc, Variant::TrailingHole);
  BOOST_REQUIRE_EQUAL(tracklets.unsmoothed.size(), 1u);
  BOOST_CHECK_EQUAL(tracklets.unsmoothed.front(), tracklets.first.tipIndex());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(JoinSmoothsStemStatesBeforeShared, smoother_t,
                              Smoothers) {
  // a hole before the first measurement gets a smoothing step from the
  // corrected shared state and matches the joint posterior; the link then
  // takes it off the chain, because the second tracklet covers that side
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto tracklets = joinAndCheck<smoother_t>(tc, Variant::StemHole);
  const auto nodes = makeModel(Variant::StemHole);
  const auto [mean, cov] = densePosterior(nodes);
  std::size_t stemIndex = kTrackIndexInvalid;
  for (const auto& [index, node] : tracklets.nodeOfState) {
    if (node == 7) {
      stemIndex = index;
    }
  }
  BOOST_REQUIRE_NE(stemIndex, kTrackIndexInvalid);
  checkStateAgainstDense(tc.trackStateContainer().getTrackState(stemIndex), 7,
                         mean, cov, 1e-9);
  // off the chain: the shared state now links into the second tracklet
  BOOST_CHECK_NE(tracklets.shared.previous(), stemIndex);
}

BOOST_AUTO_TEST_CASE(JoinCopyLeavesInputsUntouched) {
  const auto nodes = makeModel(Variant::Base);
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto tracklets = buildTracklets(tc, nodes, Variant::Base);
  const BoundVector sharedBefore = tracklets.shared.smoothed();
  const std::size_t firstStates = tracklets.first.nTrackStates();
  const std::size_t secondStates = tracklets.second.nTrackStates();
  const std::size_t tracksBefore = tc.size();

  TrackletJoiner joiner;
  auto joined = joiner.joinCopy(gctx, tc, tracklets.first, tracklets.shared,
                                tracklets.second);
  BOOST_REQUIRE(joined.ok());

  // the copy has new state indices: the joined chain reads the first
  // tracklet from its tip, then the second tracklet from its start
  std::vector<std::size_t> order;
  for (const auto& ts : tracklets.first.trackStatesReversed()) {
    order.push_back(tracklets.nodeOfState.at(ts.index()));
  }
  std::vector<std::size_t> secondOrder;
  for (const auto& ts : tracklets.second.trackStatesReversed()) {
    secondOrder.push_back(tracklets.nodeOfState.at(ts.index()));
  }
  order.insert(order.end(), secondOrder.rbegin(), secondOrder.rend());
  std::map<std::size_t, std::size_t> nodeOfCopy;
  std::size_t pos = 0;
  for (const auto& ts : joined->trackStatesReversed()) {
    nodeOfCopy[ts.index()] = order.at(pos++);
  }
  checkJoinedAgainstDense(*joined, nodeOfCopy, {}, nodes, 1e-9, nodes.size());

  CHECK_CLOSE_ABS(tracklets.shared.smoothed(), sharedBefore, 0.);
  BOOST_CHECK_EQUAL(tracklets.first.nTrackStates(), firstStates);
  BOOST_CHECK_EQUAL(tracklets.second.nTrackStates(), secondStates);
  // one joined track added, the shell of the second copy removed
  BOOST_CHECK_EQUAL(tc.size(), tracksBefore + 1);
  BOOST_CHECK_EQUAL(joined->nTrackStates(), firstStates + secondStates);
}

BOOST_AUTO_TEST_CASE(ConditioningRejectsStateOffTheChain) {
  const auto nodes = makeModel(Variant::Base);
  Container tc{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto tracklets = buildTracklets(tc, nodes, Variant::Base);
  auto res = conditionSmoothedStatesOnSharedState(
      tc.trackStateContainer(), tracklets.second.tipIndex(),
      tracklets.first.tipIndex(), tracklets.shared.smoothed(),
      tracklets.shared.smoothedCovariance());
  BOOST_CHECK(!res.ok());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
