// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/KalmanVertexTrackUpdater.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {
namespace VertexSmoothing {

/// @brief Updates all tracks at vertex
/// with knowledge of the vertex position
///
/// @param gctx The Geometry Context
/// @param vtx The vertex
///
/// @tparam input_track_t Track object type
template <typename input_track_t>
static Result<void> smoothVertexSequentially(const GeometryContext& gctx,
                                             Vertex<input_track_t>* vtx) {
  if (vtx == nullptr) {
    return VertexingError::EmptyInput;
  }

  std::vector<TrackAtVertex<input_track_t>> tracks = vtx->tracks();
  for (auto& trk : tracks) {
    // update trk
    auto res = KalmanVertexTrackUpdater::update<input_track_t>(gctx, trk, vtx);
    if (!res.ok()) {
      return res.error();
    }
  }

  vtx->setTracksAtVertex(tracks);

  return {};
}

}  // Namespace VertexSmoothing
}  // Namespace Acts
