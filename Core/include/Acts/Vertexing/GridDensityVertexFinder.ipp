// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <int mainGridSize, int trkGridSize, typename vfitter_t>
auto Acts::GridDensityVertexFinder<mainGridSize, trkGridSize, vfitter_t>::find(
    const std::vector<InputTrack>& trackVector,
    const VertexingOptions& vertexingOptions, IVertexFinder::State& state) const
    -> Result<std::vector<Vertex>> {

    auto& thisState = state.as<State>();
  // Remove density contributions from tracks removed from track collection
  if (m_cfg.cacheGridStateForTrackRemoval && thisState.isInitialized &&
      !thisState.tracksToRemove.empty()) {
    // Bool to check if removable tracks, that pass selection, still exist
    bool couldRemoveTracks = false;
    for (auto trk : thisState.tracksToRemove) {
      if (!thisState.trackSelectionMap.at(trk)) {
        // Track was never added to grid, so cannot remove it
        continue;
      }
      couldRemoveTracks = true;
      auto binAndTrackGrid = thisState.binAndTrackGridMap.at(trk);
      m_cfg.gridDensity.removeTrackGridFromMainGrid(
          binAndTrackGrid.first, binAndTrackGrid.second, thisState.mainGrid);
    }
    if (!couldRemoveTracks) {
      // No tracks were removed anymore
      // Return empty seed, i.e. vertex at constraint position
      // (Note: Upstream finder should check for this break condition)
      std::vector<Vertex> seedVec{vertexingOptions.constraint};
      return seedVec;
    }
  } else {
    thisState.mainGrid = MainGridVector::Zero();
    // Fill with track densities
    for (auto trk : trackVector) {
      const BoundTrackParameters& trkParams = m_cfg.extractParameters(trk);
      // Take only tracks that fulfill selection criteria
      if (!doesPassTrackSelection(trkParams)) {
        if (m_cfg.cacheGridStateForTrackRemoval) {
          thisState.trackSelectionMap[trk] = false;
        }
        continue;
      }
      auto binAndTrackGrid =
          m_cfg.gridDensity.addTrack(trkParams, thisState.mainGrid);
      // Cache track density contribution to main grid if enabled
      if (m_cfg.cacheGridStateForTrackRemoval) {
        thisState.binAndTrackGridMap[trk] = binAndTrackGrid;
        thisState.trackSelectionMap[trk] = true;
      }
    }
    thisState.isInitialized = true;
  }

  double z = 0;
  double width = 0;
  if (thisState.mainGrid != MainGridVector::Zero()) {
    if (!m_cfg.estimateSeedWidth) {
      // Get z value of highest density bin
      auto maxZres = m_cfg.gridDensity.getMaxZPosition(thisState.mainGrid);

      if (!maxZres.ok()) {
        return maxZres.error();
      }
      z = *maxZres;
    } else {
      // Get z value of highest density bin and width
      auto maxZres = m_cfg.gridDensity.getMaxZPositionAndWidth(thisState.mainGrid);

      if (!maxZres.ok()) {
        return maxZres.error();
      }
      z = (*maxZres).first;
      width = (*maxZres).second;
    }
  }

  // Construct output vertex
  Vector3 seedPos = vertexingOptions.constraint.position() + Vector3(0., 0., z);

  Vertex returnVertex = Vertex(seedPos);

  SquareMatrix4 seedCov = vertexingOptions.constraint.fullCovariance();

  if (width != 0.) {
    // Use z-constraint from seed width
    seedCov(2, 2) = width * width;
  }

  returnVertex.setFullCovariance(seedCov);

  std::vector<Vertex> seedVec{returnVertex};

  return seedVec;
}

template <int mainGridSize, int trkGridSize, typename vfitter_t>
auto Acts::GridDensityVertexFinder<mainGridSize, trkGridSize, vfitter_t>::
    doesPassTrackSelection(const BoundTrackParameters& trk) const -> bool {
  // Get required track parameters
  const double d0 = trk.parameters()[BoundIndices::eBoundLoc0];
  const double z0 = trk.parameters()[BoundIndices::eBoundLoc1];
  // Get track covariance
  if (!trk.covariance().has_value()) {
    return false;
  }
  const auto perigeeCov = *(trk.covariance());
  const double covDD =
      perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc0);
  const double covZZ =
      perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundLoc1);
  const double covDZ =
      perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc1);
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
