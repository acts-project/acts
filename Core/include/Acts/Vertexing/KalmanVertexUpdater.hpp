// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace KalmanVertexUpdater {

/// KalmanVertexUpdater
///
/// @brief adds or removes track from or updates current vertex
/// Based on R. Frühwirth et al.
/// Vertex reconstruction and track bundling at the lep collider using
/// robust Algorithms Computer Physics Comm.: 96 (1996) 189, chapter 2.1

/// Cache object to store matrix information
struct Cache {
  Vector4 newVertexPos = Vector4::Zero();
  SquareMatrix4 newVertexCov = SquareMatrix4::Zero();
  SquareMatrix4 newVertexWeight = SquareMatrix4::Zero();
  SquareMatrix4 oldVertexWeight = SquareMatrix4::Zero();
  // W_k from the reference
  SquareMatrix3 wMat = SquareMatrix3::Zero();
};

/// @brief Updates vertex with knowledge of new track
/// @note KalmanVertexUpdater updates the vertex w.r.t. the
/// newly given track, but does NOT add the track to the
/// TrackAtVertex list. Has to be done manually after calling
/// this method
///
/// @param vtx Vertex to be updated
/// @param trk Track to be used for updating the vertex
template <typename input_track_t>
void updateVertexWithTrack(Vertex<input_track_t>& vtx,
                           TrackAtVertex<input_track_t>& trk);

/// @brief Calculates updated vertex position and covariance when
/// adding/removing linTrack and saves the results in cache
/// @param vtx Vertex
/// @param linTrack Linearized version of track to be added or removed
/// @param trackWeight Track weight
/// @param sign +1 (add track) or -1 (remove track)
/// @note Tracks are removed during the smoothing procedure to compute
/// the chi2 of the track wrt the updated vertex position
/// @param[out] cache A cache to store the results of this function
template <typename input_track_t>
void calculateUpdate(const Acts::Vertex<input_track_t>& vtx,
                     const Acts::LinearizedTrack& linTrack,
                     const double& trackWeight, int sign, Cache& cache);

namespace detail {

/// @brief Takes old and new vtx and calculates position chi2
///
/// @param oldVtx Old vertex
/// @param cache Cache containing updated vertex position
///
/// @return Chi2
template <typename input_track_t>
double vertexPositionChi2(const Vertex<input_track_t>& oldVtx,
                          const Cache& cache);

/// @brief Calculates chi2 of refitted track parameters
/// w.r.t. updated vertex
///
/// @param linTrack Linearized version of track
/// @param cache Cache containing some quantities needed in
/// this function
///
/// @return Chi2
template <typename input_track_t>
double trackParametersChi2(const LinearizedTrack& linTrack, const Cache& cache);

/// @brief Adds or removes (depending on `sign`) tracks from vertex
/// and updates the vertex
///
/// @param vtx Vertex to be updated
/// @param trk Track to be added to/removed from vtx
/// @param sign +1 (add track) or -1 (remove track)
/// @note Tracks are removed during the smoothing procedure to compute
/// the chi2 of the track wrt the updated vertex position
template <typename input_track_t>
void update(Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk,
            int sign);
}  // Namespace detail

}  // Namespace KalmanVertexUpdater
}  // Namespace Acts

#include "KalmanVertexUpdater.ipp"
