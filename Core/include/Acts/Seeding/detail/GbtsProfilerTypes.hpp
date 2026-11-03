// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>

namespace Acts::Experimental::detail {

/// Node cuts, applied when a space point is inserted.
enum class GbtsNodeCut : std::uint8_t {
  /// Pixel endcap cluster wider than `maxEndcapClusterWidth`.
  EndcapClusterWidth,
};

/// Doublet cuts, applied by `GbtsGraphBuilder` to the pairs in a phi window.
enum class GbtsDoubletCut : std::uint8_t {
  /// The outer node has no outgoing link to cut its z0 against.
  IsolatedOuterNode,
  /// The radial separation is below `minDeltaRadius`.
  MinDeltaRadius,
  /// The pair's direction does not cross one of its strips.
  StripCrossing,
  /// |cot(theta)| above `maxAbsTau`.
  MaxAbsTau,
  /// Outside the inner node's cluster width tau window.
  InnerTauWindow,
  /// Outside the outer node's cluster width tau window.
  OuterTauWindow,
  /// z0 not in the outer node's z0 histogram.
  Z0Histogram,
  /// z0 outside [`minZ0`, `maxZ0`].
  Z0Range,
  /// z at `maxOuterRadius` outside the RoI.
  ZAtOuterRadius,
  /// |curvature| above the curvature cut.
  Curvature,
  /// No incoming edge of the outer node matched in tau.
  MatchBeforeCreate,
  /// The edge storage was full.
  EdgeLimit,
};

/// Link cuts, applied by `GbtsGraphBuilder` to a new edge and an incoming
/// edge of its outer node.
enum class GbtsTripletCut : std::uint8_t {
  /// |tau ratio| above the looser of the two layer pairs' tau ratio cuts.
  TauRatio,
  /// |dPhi| above `cutDPhiMax`.
  DeltaPhi,
  /// |dCurvature| above `cutDCurvMax`.
  DeltaCurvature,
  /// The tau ratio, dPhi and dCurvature, each divided by its cut, add up to
  /// more than `sumRatioCutMax`.
  SumRatio,
  /// The barrel triplet validation could not fit the three points.
  ValidationDegenerate,
  /// The barrel triplet validation found |d0| above `d0Max`.
  ValidationD0,
  /// The barrel triplet validation found pT below its minimum.
  ValidationPt,
  /// The barrel triplet validation found a high pT triplet with a loose tau
  /// ratio.
  ValidationHighPtTauRatio,
  /// Passed every cut, but the outer edge already has
  /// `kGbtsMaxEdgeNeighbours` links.
  NeighbourLimit,
};

/// Tracking filter cuts, applied while a path is extended by one edge, and
/// the post-fit cut on the complete chain.
enum class GbtsFilterCut : std::uint8_t {
  /// chi2 in the transverse plane above `maxDChi2X`.
  Chi2X,
  /// chi2 in r-z above `maxDChi2Y`.
  Chi2Y,
  /// |curvature| of the updated state above `maxCurvature`.
  Curvature,
  /// |z0| of the updated state above `maxZ0`.
  Z0,
  /// |curvature| of a complete chain times `postFitCurvatureFactor` above
  /// `maxCurvature`.
  PostFitCurvature,
};

/// What became of an edge as the head of a path.
enum class GbtsHeadOutcome : std::uint8_t {
  /// Its best path is a seed candidate.
  Candidate,
  /// It lies in a chain longer than `maxChainLength`.
  ChainTooLong,
  /// The longest path from it is shorter than the minimum seed length.
  ChainTooShort,
  /// The tracking filter rejected the edge itself.
  FilterRejected,
  /// The tracking filter cut every path from it short of the minimum, or the
  /// post-fit curvature cut rejected those that got there.
  NoPathSurvived,
};

/// What became of a seed candidate.
enum class GbtsCandidateFate : std::uint8_t {
  /// Output as a seed.
  Seed,
  /// Output as two seeds, each dropping one space point.
  Split,
  /// Removed as a clone of a better candidate.
  Clone,
  /// Dropped for having fewer than three space points.
  TooFewSpacePoints,
};

/// Why a check rejected a candidate, filled in by the check for the caller
/// to record.
/// @tparam Cut The cut enumeration of the check
template <typename Cut>
struct GbtsCutFailure {
  /// The cut it failed.
  Cut cut{};
  /// The value cut on.
  float value{};
  /// The threshold it failed, zero for a cut without one.
  float threshold{};
};
/// Why the barrel triplet validation rejected a triplet.
using GbtsTripletFailure = GbtsCutFailure<GbtsTripletCut>;
/// Why the tracking filter rejected an edge.
using GbtsFilterFailure = GbtsCutFailure<GbtsFilterCut>;

}  // namespace Acts::Experimental::detail
