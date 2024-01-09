// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace KalmanVertexTrackUpdater {

/// KalmanVertexTrackUpdater
/// Based on
/// Ref. (1):
/// R. Frühwirth et al.
/// Vertex reconstruction and track bundling at the lep collider using robust
/// algorithms
/// Computer Physics Comm.: 96 (1996) 189
/// Chapter 2.1

/// @brief Refits a single track with the knowledge of
/// the vertex it has originated from
///
/// @param track Track to update
/// @param vtxPosFull full vertex position
/// @param vtxCovFull full vertex covariance matrix
/// @param nDimVertex number of dimensions of the vertex
void update(TrackAtVertexRef track, const Vector4& vtxPosFull,
            const SquareMatrix4& vtxCovFull, unsigned int nDimVertex);

namespace detail {

/// @brief Calculates a covariance matrix for the refitted track parameters
///
/// @param wMat W_k matrix from Ref. (1)
/// @param crossCovVP Cross-covariance matrix between vertex position and track
/// momentum
/// @param vtxCov Vertex covariance matrix
/// @param newTrkParams Refitted track parameters
template <unsigned int nDimVertex>
inline BoundMatrix calculateTrackCovariance(
    const SquareMatrix3& wMat, const ActsMatrix<nDimVertex, 3>& crossCovVP,
    const ActsSquareMatrix<nDimVertex>& vtxCov,
    const BoundVector& newTrkParams);

}  // Namespace detail

}  // Namespace KalmanVertexTrackUpdater
}  // Namespace Acts
