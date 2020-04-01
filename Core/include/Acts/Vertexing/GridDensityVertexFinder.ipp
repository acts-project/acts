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
    const VertexingOptions<InputTrack_t>& vertexingOptions) const
    -> Result<std::vector<Vertex<InputTrack_t>>> {
  // Use this method only if no caching of density values is desired
  assert(m_cfg.cacheGridStateForTrackRemoval == false);

  // Dummy state
  State dummyState;
  return find(trackVector, vertexingOptions, dummyState);
}

template <int mainGridSize, int trkGridSize, typename vfitter_t>
auto Acts::GridDensityVertexFinder<mainGridSize, trkGridSize, vfitter_t>::find(
    const std::vector<const InputTrack_t*>& trackVector,
    const VertexingOptions<InputTrack_t>& vertexingOptions, State& state) const
    -> Result<std::vector<Vertex<InputTrack_t>>> {
  // Remove density contributions from tracks removed from track collection
  if (m_cfg.cacheGridStateForTrackRemoval && state.isInitialized &&
      !state.tracksToRemove.empty()) {
    for (const auto& trk : state.tracksToRemove) {
      auto binAndTrackGrid = state.binAndTrackGridMap.at(trk);
      m_cfg.gridDensity.removeTrackGridFromMainGrid(
          binAndTrackGrid.first, binAndTrackGrid.second, state.mainGrid);
    }
  } else {
    // Fill with track densities
    for (const auto& trk : trackVector) {
      auto binAndTrackGrid =
          m_cfg.gridDensity.addTrack(m_extractParameters(*trk), state.mainGrid);
      // Cache track density contribution to main grid if enabled
      if (m_cfg.cacheGridStateForTrackRemoval) {
        state.binAndTrackGridMap[trk] = binAndTrackGrid;
      }
    }
    state.isInitialized = true;
  }

  // Get z value of highest density bin
  auto maxZres = m_cfg.gridDensity.getMaxZPosition(state.mainGrid);

  if (!maxZres.ok()) {
    return maxZres.error();
  }

  // Construct output vertex
  Vector3D seedPos = Vector3D(0., 0., *maxZres);

  Vertex<InputTrack_t> returnVertex = Vertex<InputTrack_t>(seedPos);

  ActsSymMatrixD<3> seedCov = vertexingOptions.vertexConstraint.covariance();

  returnVertex.setCovariance(seedCov);

  std::vector<Vertex<InputTrack_t>> seedVec{returnVertex};

  return seedVec;
}
