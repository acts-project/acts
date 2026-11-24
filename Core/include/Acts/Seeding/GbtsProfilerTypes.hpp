// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <string_view>

namespace Acts::Experimental {

/// Node cuts, applied when a space point is inserted.
namespace GbtsNodeCut {
/// Pixel endcap cluster wider than `maxEndcapClusterWidth`.
inline constexpr std::string_view EndcapClusterWidth = "EndcapClusterWidth";

/// Every node cut.
inline constexpr std::array all = {EndcapClusterWidth};
}  // namespace GbtsNodeCut

/// Doublet cuts, applied by `GbtsGraphBuilder` to the pairs in a phi window.
namespace GbtsDoubletCut {
/// The outer node has no outgoing link to cut its z0 against.
inline constexpr std::string_view IsolatedOuterNode = "IsolatedOuterNode";
/// With `dropDeadEnds`, the outer node has no edges for the edge to link to.
inline constexpr std::string_view DeadEndOuterNode = "DeadEndOuterNode";
/// With `pruneShortChains`, the outer node's edges are on chains too short
/// for an edge to them to reach a seed: the value is the level such an edge
/// would have, the threshold the level its bin needs.
inline constexpr std::string_view ShortChainOuterNode = "ShortChainOuterNode";
/// The radial separation is below `minDeltaRadius`.
inline constexpr std::string_view MinDeltaRadius = "MinDeltaRadius";
/// The pair's direction does not cross one of its strips.
inline constexpr std::string_view StripCrossing = "StripCrossing";
/// |cot(theta)| above `maxAbsTau`.
inline constexpr std::string_view MaxAbsTau = "MaxAbsTau";
/// Outside the inner node's cluster width tau window.
inline constexpr std::string_view InnerTauWindow = "InnerTauWindow";
/// Outside the outer node's cluster width tau window.
inline constexpr std::string_view OuterTauWindow = "OuterTauWindow";
/// z0 not in the outer node's z0 histogram.
inline constexpr std::string_view Z0Histogram = "Z0Histogram";
/// z0 outside the geometry's z0 window.
inline constexpr std::string_view Z0Range = "Z0Range";
/// z at `maxOuterRadius` outside the RoI.
inline constexpr std::string_view ZAtOuterRadius = "ZAtOuterRadius";
/// |curvature| above the curvature cut.
inline constexpr std::string_view Curvature = "Curvature";
/// No incoming edge of the outer node matched in tau.
inline constexpr std::string_view MatchBeforeCreate = "MatchBeforeCreate";
/// The edge storage was full.
inline constexpr std::string_view EdgeLimit = "EdgeLimit";
/// With `dropDeadEnds`, the edge linked to no outer edge.
inline constexpr std::string_view DeadEnd = "DeadEnd";
/// With `pruneShortChains`, the edge's level, the longest chain outwards,
/// is below the level its bin needs.
inline constexpr std::string_view ShortChain = "ShortChain";

/// Every doublet cut.
inline constexpr std::array all = {
    IsolatedOuterNode, DeadEndOuterNode, ShortChainOuterNode,
    MinDeltaRadius,    StripCrossing,    MaxAbsTau,
    InnerTauWindow,    OuterTauWindow,   Z0Histogram,
    Z0Range,           ZAtOuterRadius,   Curvature,
    MatchBeforeCreate, EdgeLimit,        DeadEnd,
    ShortChain};
}  // namespace GbtsDoubletCut

/// Link cuts, applied by `GbtsGraphBuilder` to a new edge and an incoming
/// edge of its outer node.
namespace GbtsTripletCut {
/// |tau ratio| above the looser of the two layer pairs' tau ratio cuts.
inline constexpr std::string_view TauRatio = "TauRatio";
/// |dPhi| above `cutDPhiMax`.
inline constexpr std::string_view DeltaPhi = "DeltaPhi";
/// |dCurvature| above `cutDCurvMax`.
inline constexpr std::string_view DeltaCurvature = "DeltaCurvature";
/// The tau ratio, dPhi and dCurvature, each divided by its cut, add up to
/// more than `sumRatioCutMax`.
inline constexpr std::string_view SumRatio = "SumRatio";
/// The barrel triplet validation could not fit the three points.
inline constexpr std::string_view ValidationDegenerate = "ValidationDegenerate";
/// The barrel triplet validation found |d0| above `d0Max`.
inline constexpr std::string_view ValidationD0 = "ValidationD0";
/// The barrel triplet validation found pT below its minimum.
inline constexpr std::string_view ValidationPt = "ValidationPt";
/// The barrel triplet validation found a high pT triplet with a loose tau
/// ratio.
inline constexpr std::string_view ValidationHighPtTauRatio =
    "ValidationHighPtTauRatio";

/// Every triplet cut.
inline constexpr std::array all = {TauRatio,
                                   DeltaPhi,
                                   DeltaCurvature,
                                   SumRatio,
                                   ValidationDegenerate,
                                   ValidationD0,
                                   ValidationPt,
                                   ValidationHighPtTauRatio};
}  // namespace GbtsTripletCut

/// Tracking filter cuts, applied while a path is extended by one edge, and
/// the post-fit cut on the complete chain.
namespace GbtsFilterCut {
/// chi2 in the transverse plane above `maxDChi2X`.
inline constexpr std::string_view Chi2X = "Chi2X";
/// chi2 in r-z above `maxDChi2Y`.
inline constexpr std::string_view Chi2Y = "Chi2Y";
/// |curvature| of the updated state above `maxCurvature`.
inline constexpr std::string_view Curvature = "Curvature";
/// z0 of the updated state outside the geometry's z0 window.
inline constexpr std::string_view Z0 = "Z0";
/// |curvature| of a complete chain times `postFitCurvatureFactor` above
/// `maxCurvature`.
inline constexpr std::string_view PostFitCurvature = "PostFitCurvature";

/// Every tracking filter cut.
inline constexpr std::array all = {Chi2X, Chi2Y, Curvature, Z0,
                                   PostFitCurvature};
}  // namespace GbtsFilterCut

/// What became of an edge as the head of a path.
namespace GbtsHeadOutcome {
/// Its best path is a seed candidate.
inline constexpr std::string_view Candidate = "Candidate";
/// It lies in a chain longer than `maxChainLength`.
inline constexpr std::string_view ChainTooLong = "ChainTooLong";
/// The longest path from it is shorter than the minimum seed length, or it
/// has no inner neighbour, so no path from it reaches three space points.
inline constexpr std::string_view ChainTooShort = "ChainTooShort";
/// The tracking filter rejected the edge itself.
inline constexpr std::string_view FilterRejected = "FilterRejected";
/// The tracking filter cut every path from it short of the minimum, or the
/// post-fit curvature cut rejected those that got there.
inline constexpr std::string_view NoPathSurvived = "NoPathSurvived";

/// Every path head outcome.
inline constexpr std::array all = {Candidate, ChainTooLong, ChainTooShort,
                                   FilterRejected, NoPathSurvived};
}  // namespace GbtsHeadOutcome

/// What became of a seed candidate.
namespace GbtsCandidateFate {
/// Output as a seed.
inline constexpr std::string_view Seed = "Seed";
/// Output as two seeds, each dropping one space point.
inline constexpr std::string_view Split = "Split";
/// Removed as a clone of a better candidate.
inline constexpr std::string_view Clone = "Clone";
/// Dropped for having fewer than three space points.
inline constexpr std::string_view TooFewSpacePoints = "TooFewSpacePoints";

/// Every candidate fate.
inline constexpr std::array all = {Seed, Split, Clone, TooFewSpacePoints};
}  // namespace GbtsCandidateFate

}  // namespace Acts::Experimental
