// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
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
/// @tparam nDimVertex number of dimensions of the vertex, i.e., be 4
/// when using time information and 3 otherwise
template <unsigned int nDimVertex>
struct MatrixCache {
  ActsVector<nDimVertex> newVertexPos = ActsVector<nDimVertex>::Zero();
  ActsSquareMatrix<nDimVertex> newVertexCov =
      ActsSquareMatrix<nDimVertex>::Zero();
  ActsSquareMatrix<nDimVertex> newVertexWeight =
      ActsSquareMatrix<nDimVertex>::Zero();
  ActsSquareMatrix<nDimVertex> oldVertexWeight =
      ActsSquareMatrix<nDimVertex>::Zero();
  SquareMatrix3 momWeightInv = SquareMatrix3::Zero();
};

/// @brief Updates vertex with knowledge of new track
/// @note KalmanVertexUpdater updates the vertex w.r.t. the
/// newly given track, but does NOT add the track to the
/// TrackAtVertex list. Has to be done manually after calling
/// this method
///
/// @tparam input_track_t track parameter type
/// @tparam nDimVertex number of dimensions of the vertex, i.e., be 4
/// when using time information and 3 otherwise
/// @param vtx Vertex to be updated
/// @param trk Track to be used for updating the vertex
/// @param useTimeInformation boolean indicating whether we
/// use time information in the vertex fit
template <typename input_track_t, unsigned int nDimVertex>
void updateVertexWithTrack(Vertex<input_track_t>& vtx,
                           TrackAtVertex<input_track_t>& trk);

/// @brief Updates vertex position
///
/// @tparam input_track_t track parameter type
/// @tparam nDimVertex number of dimensions of the vertex, i.e., be 4
/// when using time information and 3 otherwise
/// @param vtx Old vertex
/// @param linTrack Linearized version of track to be added or removed
/// @param trackWeight Track weight
/// @param sign +1 (add track) or -1 (remove track)
/// @param[out] matrixCache A cache to store matrix information
template <typename input_track_t, unsigned int nDimVertex>
void updatePosition(const Acts::Vertex<input_track_t>& vtx,
                    const Acts::LinearizedTrack& linTrack, double trackWeight,
                    int sign, MatrixCache<nDimVertex>& matrixCache);

namespace detail {

/// @brief Takes old and new vtx and calculates position chi2
///
/// @tparam input_track_t track parameter type
/// @tparam nDimVertex number of dimensions of the vertex, i.e., be 4
/// when using time information and 3 otherwise
/// @param oldVtx Old vertex
/// @param matrixCache A cache to store matrix information
///
/// @return Chi2
template <typename input_track_t, unsigned int nDimVertex>
double vertexPositionChi2(const Vertex<input_track_t>& oldVtx,
                          const MatrixCache<nDimVertex>& matrixCache);

/// @brief Calculates chi2 of refitted track parameters
/// w.r.t. updated vertex
///
/// @tparam input_track_t track parameter type
/// @tparam nDimVertex number of dimensions of the vertex, i.e., be 4
/// when using time information and 3 otherwise
/// @param linTrack Linearized version of track
/// @param matrixCache A cache to store matrix information
///
/// @return Chi2
template <typename input_track_t, unsigned int nDimVertex>
double trackParametersChi2(const LinearizedTrack& linTrack,
                           const MatrixCache<nDimVertex>& matrixCache);

/// @brief Adds or removes (depending on `sign`) tracks from vertex
/// and updates the vertex
///
/// @tparam input_track_t track parameter type
/// @tparam nDimVertex number of dimensions of the vertex, i.e., be 4
/// when using time information and 3 otherwise
/// @param vtx Vertex to be updated
/// @param trk Track to be added to/removed from vtx
/// @param sign +1 (add track) or -1 (remove track)
template <typename input_track_t, unsigned int nDimVertex>
void update(Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk,
            int sign);
}  // Namespace detail

}  // Namespace KalmanVertexUpdater
}  // Namespace Acts

#include "KalmanVertexUpdater.ipp"
