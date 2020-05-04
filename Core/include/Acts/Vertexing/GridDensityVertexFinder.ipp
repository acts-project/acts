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
      if (not state.trackSelectionMap.at(trk)) {
        // Track was never added to grid, so cannot remove it
        continue;
      }
      auto binAndTrackGrid = state.binAndTrackGridMap.at(trk);
      m_cfg.gridDensity.removeTrackGridFromMainGrid(
          binAndTrackGrid.first, binAndTrackGrid.second, state.mainGrid);
    }
  } else {
    state.mainGrid = ActsVectorF<mainGridSize>::Zero();
    // Fill with track densities
    for (auto trk : trackVector) {
      const BoundParameters& trkParams = m_extractParameters(*trk);
      // Take only tracks that fulfill selection criteria
      if (not doesPassTrackSelection(trkParams)) {
        if (m_cfg.cacheGridStateForTrackRemoval) {
          state.trackSelectionMap[trk] = false;
        }
        continue;
      }
      auto binAndTrackGrid =
          m_cfg.gridDensity.addTrack(trkParams, state.mainGrid);
      // Cache track density contribution to main grid if enabled
      if (m_cfg.cacheGridStateForTrackRemoval) {
        state.binAndTrackGridMap[trk] = binAndTrackGrid;
        state.trackSelectionMap[trk] = true;
      }
    }
    state.isInitialized = true;
  }

  double z = 0;
  double width = 0;
  if (state.mainGrid != ActsVectorF<mainGridSize>::Zero()) {
    if (not m_cfg.estimateSeedWidth) {
      // Get z value of highest density bin
      auto maxZres = m_cfg.gridDensity.getMaxZPosition(state.mainGrid);

      if (!maxZres.ok()) {
        return maxZres.error();
      }
      z = *maxZres;
    } else {
      // Get z value of highest density bin and width
      auto maxZres = m_cfg.gridDensity.getMaxZPositionAndWidth(state.mainGrid);

      if (!maxZres.ok()) {
        return maxZres.error();
      }
      z = (*maxZres).first;
      width = (*maxZres).second;
    }
  }

  // Construct output vertex
  Vector3D seedPos =
      vertexingOptions.vertexConstraint.position() + Vector3D(0., 0., z);

  Vertex<InputTrack_t> returnVertex = Vertex<InputTrack_t>(seedPos);

  ActsSymMatrixD<4> seedCov =
      vertexingOptions.vertexConstraint.fullCovariance();

  if (width != 0.) {
    // Use z-constraint from seed width
    seedCov(2, 2) = width * width;
  }

  returnVertex.setFullCovariance(seedCov);

  std::vector<Vertex<InputTrack_t>> seedVec{returnVertex};

  return seedVec;
}

template <int mainGridSize, int trkGridSize, typename vfitter_t>
auto Acts::GridDensityVertexFinder<mainGridSize, trkGridSize, vfitter_t>::
    doesPassTrackSelection(const BoundParameters& trk) const -> bool {
  // Get required track parameters
  const double d0 = trk.parameters()[ParID_t::eLOC_D0];
  const double z0 = trk.parameters()[ParID_t::eLOC_Z0];
  // Get track covariance
  const auto perigeeCov = *(trk.covariance());
  const double covDD = perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_D0);
  const double covZZ = perigeeCov(ParID_t::eLOC_Z0, ParID_t::eLOC_Z0);
  const double covDZ = perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_Z0);
  const double covDeterminant = covDD * covZZ - covDZ * covDZ;

  // Do track selection based on track cov matrix and d0SignificanceCut
  if ((covDD <= 0) || (d0 * d0 / covDD > m_cfg.d0SignificanceCut) ||
      (covZZ <= 0) || (covDeterminant <= 0)) {
    return false;
  }

  // Calculate track density quantities
  double constantTerm =
      -(d0 * d0 * covZZ + z0 * z0 * covDD + 2. * d0 * z0 * covDZ) /
      (2. * covDeterminant);
  const double linearTerm = (d0 * covDZ + z0 * covDD) / covDeterminant;
  const double quadraticTerm = -covDD / (2. * covDeterminant);
  double discriminant =
      linearTerm * linearTerm -
      4. * quadraticTerm * (constantTerm + 2. * m_cfg.z0SignificanceCut);
  if (discriminant < 0) {
    return false;
  }
  return true;
}
