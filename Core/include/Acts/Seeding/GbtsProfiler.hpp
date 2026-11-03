// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/detail/GbtsProfilerTypes.hpp"

#include <cstdint>
#include <mutex>
#include <optional>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <vector>

namespace Acts::Experimental {

class GraphBasedTrackSeeder;

/// Records what the GBTS cuts let through and what they reject, for studies of
/// where seeds are lost.
///
/// Handed to `GraphBasedTrackSeeder` at construction. Every `createSeeds` call
/// then adds one `Event`, under the number the calling thread last gave
/// `setEventNumber`, which the caller takes back out by that number or all at
/// once. A space point is always identified by the index the caller
/// inserted it with, so the records can be matched to truth afterwards.
///
/// Each candidate is recorded once, at the first cut it fails, in the order the
/// cuts are applied. A pair or triplet that is never formed, because it falls
/// outside the phi window or the layer connections, does not appear at all.
///
/// What a record failed, or became of it, is a label to read with `name`. The
/// labels' types are not part of the interface, nor are their values stable.
///
/// One busy event records on the order of ten million rejections, so a caller
/// processing many events should take each event out once it is done with it.
class GbtsProfiler {
 public:
  /// A rejected space point.
  struct Node {
    /// The space point.
    SpacePointIndex spacePoint{};
    /// The cut it failed.
    detail::GbtsNodeCut cut{};
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
    detail::GbtsDoubletCut cut{};
    /// The value cut on, e.g. the curvature for the `Curvature` cut.
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
    detail::GbtsTripletCut cut{};
    /// The value cut on, e.g. |tau ratio| for the `TauRatio` cut.
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
    detail::GbtsHeadOutcome outcome{};
    /// Edges on the longest path from it inwards, itself included.
    std::uint32_t maxPathLength{};
  };

  /// A path the tracking filter cut short, or a complete chain the post-fit
  /// curvature cut rejected, whose innermost edge then counts as the failed
  /// one.
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
    detail::GbtsFilterCut cut{};
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
    detail::GbtsCandidateFate fate{};
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
    /// The number the thread that added it last gave `setEventNumber`, zero if
    /// it gave none.
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

  /// Set the event number the calling thread's events are recorded under, until
  /// it sets another. Call it before each `createSeeds` call, on the same
  /// thread. Each thread keeps its own number, so that threads processing
  /// different events can share a profiler. Thread safe.
  /// @param eventNumber The event number
  void setEventNumber(std::uint64_t eventNumber);

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
  static std::string_view name(detail::GbtsNodeCut cut);
  /// @param cut The cut
  /// @return Its name
  static std::string_view name(detail::GbtsDoubletCut cut);
  /// @param cut The cut
  /// @return Its name
  static std::string_view name(detail::GbtsTripletCut cut);
  /// @param cut The cut
  /// @return Its name
  static std::string_view name(detail::GbtsFilterCut cut);
  /// @param outcome The outcome
  /// @return Its name
  static std::string_view name(detail::GbtsHeadOutcome outcome);
  /// @param fate The fate
  /// @return Its name
  static std::string_view name(detail::GbtsCandidateFate fate);

  /// @return The name of every node cut
  static std::vector<std::string_view> nodeCutNames();
  /// @return The name of every doublet cut
  static std::vector<std::string_view> doubletCutNames();
  /// @return The name of every triplet cut
  static std::vector<std::string_view> tripletCutNames();
  /// @return The name of every tracking filter cut
  static std::vector<std::string_view> filterCutNames();
  /// @return The name of every path head outcome
  static std::vector<std::string_view> headOutcomeNames();
  /// @return The name of every candidate fate
  static std::vector<std::string_view> candidateFateNames();

 private:
  friend class GraphBasedTrackSeeder;

  /// Add an event under the calling thread's event number. Thread safe.
  /// @param event The event to add
  void add(Event&& event);

  std::mutex m_mutex;
  std::vector<Event> m_events;
  std::unordered_map<std::thread::id, std::uint64_t> m_eventNumbers;
};

}  // namespace Acts::Experimental
