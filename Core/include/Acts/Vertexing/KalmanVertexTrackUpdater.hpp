// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
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
/// @tparam input_track_t track parameter type
/// @param track Track to update
/// @param vtx Vertex `track` belongs to
template <typename input_track_t>
void update(TrackAtVertex<input_track_t>& track,
            const Vertex<input_track_t>& vtx);

namespace detail {

/// @brief reates a new covariance matrix for the
/// refitted track parameters
///
/// @param sMat Track ovariance in momentum space
/// @param crossCovVP Cross variance between the vertex and
/// the fitted track momentum
/// @param vtxWeight Vertex weight matrix
/// @param vtxCov Vertex covariance matrix
/// @param newTrkParams New track parameter
BoundMatrix calculateTrackCovariance(const SquareMatrix3& sMat,
                                     const ActsMatrix<4, 3>& crossCovVP,
                                     const SquareMatrix4& vtxWeight,
                                     const SquareMatrix4& vtxCov,
                                     const BoundVector& newTrkParams);

}  // Namespace detail

}  // Namespace KalmanVertexTrackUpdater
}  // Namespace Acts

#include "Acts/Vertexing/KalmanVertexTrackUpdater.ipp"
