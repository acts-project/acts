// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"

#include <cstdint>
#include <mutex>
#include <optional>
#include <string_view>
#include <vector>

namespace Acts::Experimental {

/// Records what the GBTS cuts let through and what they reject, for studies of
/// where seeds are lost.
///
/// Handed to `GraphBasedTrackSeeder` at construction. Every `createSeeds` call
/// then adds one `Event`, which the caller takes back out by event number or
/// all at once. A space point is always identified by the index the caller
/// inserted it with, so the records can be matched to truth afterwards.
///
/// Each candidate is recorded once, at the first cut it fails, in the order the
/// cuts are applied. A pair or triplet that is never formed, because it falls
/// outside the phi window or the layer connections, does not appear at all.
///
/// One busy event records on the order of ten million rejections, so a caller
/// processing many events should take each event out once it is done with it.
class GbtsProfiler {
 public:
  /// Node cuts, applied when a space point is inserted.
  enum class NodeCut : std::uint8_t {
    /// Pixel endcap cluster wider than `maxEndcapClusterWidth`.
    EndcapClusterWidth,
  };

  /// Doublet cuts, applied by `GbtsGraphBuilder` to the pairs in a phi window.
  enum class DoubletCut : std::uint8_t {
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
  enum class TripletCut : std::uint8_t {
    /// |tau ratio| above the looser of the two layer pairs' tau ratio cuts.
    TauRatio,
    /// |dPhi| above `cutDPhiMax`.
    DeltaPhi,
    /// |dCurvature| above `cutDCurvMax`.
    DeltaCurvature,
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

  /// Tracking filter cuts, applied while a path is extended by one edge.
  enum class FilterCut : std::uint8_t {
    /// chi2 in the transverse plane above `maxDChi2X`.
    Chi2X,
    /// chi2 in r-z above `maxDChi2Y`.
    Chi2Y,
    /// |curvature| of the updated state above `maxCurvature`.
    Curvature,
    /// |z0| of the updated state above `maxZ0`.
    Z0,
  };

  /// What became of an edge as the head of a path.
  enum class HeadOutcome : std::uint8_t {
    /// Its best path is a seed candidate.
    Candidate,
    /// It lies in a chain longer than `maxChainLength`.
    ChainTooLong,
    /// The longest path from it is shorter than the minimum seed length.
    ChainTooShort,
    /// The tracking filter rejected the edge itself.
    FilterRejected,
    /// The tracking filter cut every path from it short of the minimum.
    NoPathSurvived,
  };

  /// What became of a seed candidate.
  enum class CandidateFate : std::uint8_t {
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
  struct Failure {
    /// The cut it failed.
    Cut cut{};
    /// The value cut on.
    float value{};
    /// The threshold it failed, zero for a cut without one.
    float threshold{};
  };
  /// Why the barrel triplet validation rejected a triplet.
  using TripletFailure = Failure<TripletCut>;
  /// Why the tracking filter rejected an edge.
  using FilterFailure = Failure<FilterCut>;

  /// A rejected space point.
  struct Node {
    /// The space point.
    SpacePointIndex spacePoint{};
    /// The cut it failed.
    NodeCut cut{};
    /// The value cut on, e.g. the cluster width.
    float value{};
    /// The threshold it failed.
    float threshold{};
  };

  /// An accepted doublet, i.e. a graph edge.
  struct Edge {
    /// The inner space point.
    SpacePointIndex inner{};
    /// The outer space point.
    SpacePointIndex outer{};
  };

  /// A rejected doublet.
  struct Doublet {
    /// The inner space point.
    SpacePointIndex inner{};
    /// The outer space point.
    SpacePointIndex outer{};
    /// The cut it failed.
    DoubletCut cut{};
    /// The value cut on, e.g. the curvature for `DoubletCut::Curvature`.
    float value{};
    /// The threshold it failed, zero for a cut without one.
    float threshold{};
  };

  /// An accepted triplet, i.e. a link between two edges.
  struct Link {
    /// The inner space point.
    SpacePointIndex inner{};
    /// The middle space point, shared by the two edges.
    SpacePointIndex middle{};
    /// The outer space point.
    SpacePointIndex outer{};
  };

  /// A rejected triplet.
  struct Triplet {
    /// The inner space point.
    SpacePointIndex inner{};
    /// The middle space point, shared by the two edges.
    SpacePointIndex middle{};
    /// The outer space point.
    SpacePointIndex outer{};
    /// The cut it failed.
    TripletCut cut{};
    /// The value cut on, e.g. |tau ratio| for `TripletCut::TauRatio`.
    float value{};
    /// The threshold it failed, zero for a cut without one.
    float threshold{};
  };

  /// An edge as the head of a path, the outermost edge of its seed candidate.
  struct Head {
    /// The inner space point of the edge.
    SpacePointIndex inner{};
    /// The outer space point of the edge.
    SpacePointIndex outer{};
    /// What became of it.
    HeadOutcome outcome{};
    /// Edges on the longest path from it inwards, itself included.
    std::uint32_t maxPathLength{};
  };

  /// A path the tracking filter cut short.
  struct FilterRejection {
    /// The inner space point of the head edge the path started from.
    SpacePointIndex headInner{};
    /// The outer space point of the head edge.
    SpacePointIndex headOuter{};
    /// The inner space point of the edge that failed to extend the path, the
    /// head edge itself if that failed.
    SpacePointIndex inner{};
    /// The outer space point of that edge.
    SpacePointIndex outer{};
    /// The cut it failed.
    FilterCut cut{};
    /// Edges on the path before the failed one.
    std::uint32_t pathLength{};
    /// The value cut on.
    float value{};
    /// The threshold it failed.
    float threshold{};
    /// The space points of the path up to and including the edge that
    /// failed, innermost first: `inner`, `outer`, ..., `headInner`,
    /// `headOuter`, `pathLength + 2` in all.
    std::vector<SpacePointIndex> spacePoints;
  };

  /// A seed candidate, i.e. the best path from a head edge.
  struct Candidate {
    /// Its space points, innermost first.
    std::vector<SpacePointIndex> spacePoints;
    /// The quality the clone removal ranks by, lower is better.
    float quality{};
    /// What became of it.
    CandidateFate fate{};
    /// Its space points that better candidates own; it is a clone when
    /// more than `shareThreshold` of `spacePoints` are.
    std::uint32_t nShared{};
    /// `hitShareThreshold`, zero for a candidate the clone removal never saw.
    float shareThreshold{};
    /// The candidate, by its index in `Event::candidates`, that owns most of
    /// its shared space points (the better one on a tie), -1 if it shares
    /// none. For a clone, the candidate it is a clone of.
    std::int32_t sharedWith{-1};
  };

  /// What one `createSeeds` call recorded.
  struct Event {
    /// Event number passed in `GraphBasedTrackSeeder::Options`.
    std::uint64_t eventNumber{};
    /// Space points the node storage rejected.
    std::vector<Node> rejectedNodes;
    /// Accepted doublets, in the order the graph builder made them.
    std::vector<Edge> edges;
    /// Rejected doublets.
    std::vector<Doublet> rejectedDoublets;
    /// Accepted triplets.
    std::vector<Link> links;
    /// Rejected triplets.
    std::vector<Triplet> rejectedTriplets;
    /// What became of each edge as a path head, in the order of `edges`.
    std::vector<Head> heads;
    /// Paths the tracking filter cut short. One edge may fail on several
    /// paths from the same head.
    std::vector<FilterRejection> filterRejections;
    /// Seed candidates, one per head with the outcome `Candidate`.
    std::vector<Candidate> candidates;
  };

  /// Add an event. Thread safe.
  /// @param event The event to add
  void add(Event&& event);

  /// Take an event out. Thread safe.
  /// @param eventNumber Event number it was recorded with
  /// @return The event, or nothing if no event has that number. Of several
  ///         with the same number, the one added first.
  std::optional<Event> take(std::uint64_t eventNumber);

  /// Take every event out. Thread safe.
  /// @return The events, in the order they were added
  std::vector<Event> takeAll();

  /// @param cut The cut
  /// @return Its name
  static std::string_view name(NodeCut cut);
  /// @param cut The cut
  /// @return Its name
  static std::string_view name(DoubletCut cut);
  /// @param cut The cut
  /// @return Its name
  static std::string_view name(TripletCut cut);
  /// @param cut The cut
  /// @return Its name
  static std::string_view name(FilterCut cut);
  /// @param outcome The outcome
  /// @return Its name
  static std::string_view name(HeadOutcome outcome);
  /// @param fate The fate
  /// @return Its name
  static std::string_view name(CandidateFate fate);

 private:
  std::mutex m_mutex;
  std::vector<Event> m_events;
};

}  // namespace Acts::Experimental
