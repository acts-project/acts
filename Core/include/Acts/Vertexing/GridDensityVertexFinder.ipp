// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <int mainGridSize, int trkGridSize, typename vfitter_t>
auto Acts::GridDensityVertexFinder<mainGridSize, trkGridSize, vfitter_t>::find(
    const std::vector<const InputTrack_t*>& trackVector,
    const VertexingOptions<InputTrack_t>& vertexingOptions, State& state) const
    -> Result<std::vector<Vertex<InputTrack_t>>> {
  // Remove density contributions from tracks removed from track collection
  if (m_cfg.cacheGridStateForTrackRemoval && state.isInitialized &&
      !state.tracksToRemove.empty()) {
    for (auto trk : state.tracksToRemove) {
      auto binAndTrackGrid = state.binAndTrackGridMap.at(trk);
      m_cfg.gridDensity.removeTrackGridFromMainGrid(
          binAndTrackGrid.first, binAndTrackGrid.second, state.mainGrid);
    }
  } else {
    state.mainGrid = ActsVectorF<mainGridSize>::Zero();
    // Fill with track densities
    for (auto trk : trackVector) {
      auto binAndTrackGrid =
          m_cfg.gridDensity.addTrack(m_extractParameters(*trk), state.mainGrid);
      // Cache track density contribution to main grid if enabled
      if (m_cfg.cacheGridStateForTrackRemoval) {
        state.binAndTrackGridMap[trk] = binAndTrackGrid;
      }
    }
    state.isInitialized = true;
  }

  double z = 0;
  if (state.mainGrid != ActsVectorF<mainGridSize>::Zero()) {
    // Get z value of highest density bin
    auto maxZres = m_cfg.gridDensity.getMaxZPosition(state.mainGrid);

    if (!maxZres.ok()) {
      return maxZres.error();
    }
    z = *maxZres;
  }

  // Construct output vertex
  Vector3D seedPos =
      vertexingOptions.vertexConstraint.position() + Vector3D(0., 0., z);

  Vertex<InputTrack_t> returnVertex = Vertex<InputTrack_t>(seedPos);

  ActsSymMatrixD<4> seedCov =
      vertexingOptions.vertexConstraint.fullCovariance();

  returnVertex.setFullCovariance(seedCov);

  std::vector<Vertex<InputTrack_t>> seedVec{returnVertex};

  return seedVec;
}
