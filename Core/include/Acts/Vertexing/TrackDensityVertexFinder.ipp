// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename vfitter_t, typename track_density_t>
auto Acts::TrackDensityVertexFinder<vfitter_t, track_density_t>::find(
    const std::vector<const InputTrack_t*>& trackVector,
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& /*state*/) const -> Result<std::vector<Vertex<InputTrack_t>>> {
  typename track_density_t::State densityState(trackVector.size());

  // Calculate z seed position
  std::pair<double, double> zAndWidth =
      m_cfg.trackDensityEstimator.globalMaximumWithWidth(
          densityState, trackVector, m_extractParameters);

  double z = zAndWidth.first;

  // Calculate seed position
  // Note: constraint position is (0,0,0) if no constraint provided
  Vector3 seedPos = vertexingOptions.constraint.position() + Vector3(0., 0., z);

  Vertex<InputTrack_t> returnVertex = Vertex<InputTrack_t>(seedPos);

  SquareMatrix3 seedCov = vertexingOptions.constraint.covariance();

  // Check if a constraint is provided and set the new z position constraint
  if (seedCov != SquareMatrix3::Zero() && std::isnormal(zAndWidth.second)) {
    seedCov(eZ, eZ) = zAndWidth.second * zAndWidth.second;
  }

  returnVertex.setCovariance(seedCov);

  std::vector<Vertex<InputTrack_t>> seedVec{returnVertex};

  return seedVec;
}
