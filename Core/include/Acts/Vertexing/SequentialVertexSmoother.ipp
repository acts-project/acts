// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t>
Acts::Result<void> Acts::SequentialVertexSmoother<input_track_t>::smooth(
    const Acts::GeometryContext& gctx, Acts::Vertex<input_track_t>* vtx) const {
  if (vtx == nullptr) {
    return VertexingError::EmptyInput;
  }

  std::vector<TrackAtVertex<input_track_t>> tracks = vtx->tracks();
  for (auto& trk : tracks) {
    // update trk
    auto res = m_cfg.trackUpdator.update(gctx, trk, vtx);
    if (!res.ok()) {
      return res.error();
    }
  }

  vtx->setTracksAtVertex(tracks);

  return {};
}
