// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>
#include <map>
#include <vector>

namespace Acts {

/// @brief Problem-level description of a single vertex candidate
///
/// This holds everything a caller has to supply about one vertex before it can
/// be fitted: the constraint to fit against, the position the seed finder
/// proposed, and the tracks assigned to it. It deliberately contains no
/// per-iteration scratch data, so it stays meaningful independently of any
/// particular fitter.
struct VertexFitCandidate {
  VertexFitCandidate() = default;

  /// Construct a candidate from a constraint and a seed position.
  /// @param constr Vertex constraint for the fitting procedure
  /// @param pos Seed position as proposed by the vertex seed finder
  VertexFitCandidate(const Acts::Vertex& constr, const Acts::Vector4& pos)
      : constraint(constr), seedPosition(pos) {}

  /// Vertex constraint for the fitting procedure
  Acts::Vertex constraint;

  /// The seed position, i.e. the first estimate of the vertex position as
  /// obtained by the vertex seed finder
  Acts::Vector4 seedPosition{Acts::Vector4::Zero()};

  /// All tracks that are currently assigned to this vertex
  std::vector<InputTrack> trackLinks;
};

/// @brief A multi-vertex fit problem
///
/// Per-event, caller-owned and fitter-agnostic description of what is to be
/// fitted: a set of vertices, the candidate description of each, and the
/// track-to-vertex association that couples them. Fitter-private scratch data
/// (linearization points, annealing state, field caches) lives in the
/// respective fitter's cache instead, so that this type can be shared between
/// fitters and inspected by callers.
struct VertexFitProblem {
  /// The vertices to be fitted
  std::vector<Vertex*> vertices;

  /// Candidate description (constraint, seed position, tracks) per vertex
  std::map<Vertex*, VertexFitCandidate> candidates;

  /// Multimap connecting tracks to all of their associated vertices
  std::multimap<InputTrack, Vertex*> trackToVertices;

  /// Track-at-vertex information for each (track, vertex) pair
  std::map<std::pair<InputTrack, Vertex*>, TrackAtVertex> tracksAtVertices;

  /// Adds a vertex to @c trackToVertices, using the tracks currently
  /// assigned to it in @c candidates
  /// @param vtx Vertex to add to the multimap along with its track associations
  void addVertexToMultiMap(Vertex& vtx) {
    for (auto trk : candidates[&vtx].trackLinks) {
      trackToVertices.emplace(trk, &vtx);
    }
  }

  /// Removes a vertex from @c trackToVertices
  /// @param vtx Vertex to remove from the multimap along with its track associations
  void removeVertexFromMultiMap(const Vertex& vtx) {
    for (auto iter = trackToVertices.begin(); iter != trackToVertices.end();) {
      if (iter->second == &vtx) {
        iter = trackToVertices.erase(iter);
      } else {
        ++iter;
      }
    }
  }

  /// Remove a vertex from the vertex collection
  /// @param vtxToRemove Vertex to remove from the collection
  /// @param logger Logger for diagnostic messages
  /// @return Result indicating success or failure of the removal operation
  Result<void> removeVertexFromCollection(Vertex& vtxToRemove,
                                          const Logger& logger) {
    auto it = std::ranges::find(vertices, &vtxToRemove);
    // Check if the value was found before erasing
    if (it == vertices.end()) {
      ACTS_ERROR("vtxToRemove is not part of the vertex collection.");
      return VertexingError::ElementNotFound;
    }
    // Erase the element if found
    vertices.erase(it);
    return {};
  }
};

}  // namespace Acts
