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
  // Construct empty 1-d grid along z-axis to be filled
  ActsVectorF<mainGridSize> mainGrid(ActsVectorF<mainGridSize>::Zero());

  // Fill with track densities
  for (const auto& trk : trackVector) {
    auto binAndTrackGrid =
        m_cfg.gridDensity.addTrack(m_extractParameters(*trk), mainGrid);
  }

  // Get z value of highest density bin
  auto maxZres = m_cfg.gridDensity.getMaxZPosition(mainGrid);

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
