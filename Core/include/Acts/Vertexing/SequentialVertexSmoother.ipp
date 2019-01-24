// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename input_track_t>
void Acts::SequentialVertexSmoother<input_track_t>::smooth(
    const Acts::GeometryContext& gctx, Acts::Vertex<input_track_t>& vtx) const {
  std::vector<TrackAtVertex<input_track_t>> tracks = vtx.tracks();
  for (auto& trk : tracks) {
    // update trk
    m_cfg.trackUpdator.update(gctx, trk, vtx);
  }

  vtx.setTracksAtVertex(tracks);
}
