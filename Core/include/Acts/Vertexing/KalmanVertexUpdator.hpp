// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

/// @class KalmanVertexUpdator
///
/// @brief adds or removes track from or updates current vertex
/// Based on R. Frühwirth et al.
/// Vertex reconstruction and track bundling at the lep collider using
/// robust Algorithms Computer Physics Comm.: 96 (1996) 189, chapter 2.1
///
/// @tparam input_track_t Track object type
template <typename input_track_t>
class KalmanVertexUpdator {
 public:
  /// Default constructor
  KalmanVertexUpdator() = default;

  /// @brief Method adding single track to vertex estimate
  /// and updates vertex
  ///
  /// @param vtx Vertex to be updated
  /// @param trk Track to be added to vtx
  void addAndUpdate(Vertex<input_track_t>* vtx,
                    TrackAtVertex<input_track_t>& trk) const;

  /// @brief Updates vertex position
  ///
  /// @param vtx Old vertex
  /// @param linTrack Linearized version of track to be added or removed
  /// @param trackWeight Track weight
  /// @param sign +1 (add track) or -1 (remove track)
  ///
  /// @return Vertex with updated position and covariance
  Vertex<input_track_t> updatePosition(const Vertex<input_track_t>* vtx,
                                       const LinearizedTrack& linTrack,
                                       double trackWeight, int sign) const;

 private:
  /// @brief Takes old and new vtx and calculates position chi2
  ///
  /// @param oldVtx Old vertex
  /// @param newVtx New vertex
  ///
  /// @return Chi2
  double vertexPositionChi2(const Vertex<input_track_t>* oldVtx,
                            const Vertex<input_track_t>* newVtx) const;

  /// @brief Calculates chi2 of refitted track parameters
  /// w.r.t. updated vertex
  ///
  /// @param vtx The already updated vertex
  /// @param linTrack Linearized version of track
  ///
  /// @return Chi2
  double trackParametersChi2(const Vertex<input_track_t>& vtx,
                             const LinearizedTrack& linTrack) const;

  /// @brief Adds or removes (depending on `sign`) tracks from vertex
  /// and updates the vertex
  ///
  /// @param vtx Vertex to be updated
  /// @param trk Track to be added to/removed from vtx
  /// @param sign +1 (add track) or -1 (remove track)
  void update(Vertex<input_track_t>* vtx, TrackAtVertex<input_track_t>& trk,
              int sign) const;

  /// @brief Removes track from vertex if it is attached
  ///
  /// @param vtx The vertex
  /// @param trk The track
  void removeTrackIf(Vertex<input_track_t>* vtx,
                     const TrackAtVertex<input_track_t>& trk) const;
};

}  // Namespace Acts

#include "KalmanVertexUpdator.ipp"