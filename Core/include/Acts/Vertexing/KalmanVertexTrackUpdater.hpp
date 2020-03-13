// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace KalmanVertexTrackUpdater {

/// KalmanVertexTrackUpdater
///
/// @brief Refits a single track with the knowledge of
/// the vertex it has originated from
///
/// @param gctx The Geometry Context
/// @param track Track to update
/// @param vtx Vertex `track` belongs to
template <typename input_track_t>
void update(const GeometryContext& gctx, TrackAtVertex<input_track_t>& track,
            const Vertex<input_track_t>& vtx);

namespace detail {

/// @brief reates a new covariance matrix for the
/// refitted track parameters
///
/// @param sMat Track ovariance in momentum space
/// @param newTrkCov New track covariance matrixs
/// @param vtxWeight Vertex weight matrix
/// @param vtxCov Vertex covariance matrix
/// @param newTrkParams New track parameter
BoundMatrix createFullTrackCovariance(
    const ActsSymMatrixD<3>& sMat,
    const ActsMatrixD<SpacePointDim, 3>& newTrkCov,
    const SpacePointSymMatrix& vtxWeight, const SpacePointSymMatrix& vtxCov,
    const BoundVector& newTrkParams);

}  // Namespace detail

}  // Namespace KalmanVertexTrackUpdater
}  // Namespace Acts

#include "Acts/Vertexing/KalmanVertexTrackUpdater.ipp"