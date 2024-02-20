// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename track_density_t>
auto Acts::TrackDensityVertexFinder<track_density_t>::find(
    const std::vector<InputTrack>& trackVector,
    const VertexingOptions& vertexingOptions,
    IVertexFinder::State& /*state*/) const -> Result<std::vector<Vertex>> {
  typename track_density_t::State densityState(trackVector.size());

  // Calculate z seed position
  std::pair<double, double> zAndWidth =
      m_cfg.trackDensityEstimator.globalMaximumWithWidth(densityState,
                                                         trackVector);

  double z = zAndWidth.first;

  // Calculate seed position
  // Note: constraint position is (0,0,0) if no constraint provided
  Vector4 seedPos =
      vertexingOptions.constraint.fullPosition() + Vector4(0., 0., z, 0.);

  Vertex returnVertex = Vertex(seedPos);

  SquareMatrix4 seedCov = vertexingOptions.constraint.fullCovariance();

  // Check if a constraint is provided and set the new z position constraint
  if (seedCov != SquareMatrix4::Zero() && std::isnormal(zAndWidth.second)) {
    seedCov(eZ, eZ) = zAndWidth.second * zAndWidth.second;
  }

  returnVertex.setFullCovariance(seedCov);

  return {returnVertex};
}
