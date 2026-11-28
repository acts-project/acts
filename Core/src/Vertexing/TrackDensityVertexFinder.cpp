// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"

#include "Acts/Vertexing/VertexingError.hpp"

Acts::Result<std::vector<Acts::Vertex>> Acts::TrackDensityVertexFinder::find(
    const std::vector<InputTrack>& trackVector,
    const VertexingOptions& vertexingOptions,
    IVertexFinder::State& /*state*/) const {
  if (vertexingOptions.fitConstraint) {
    return VertexingError::InvalidInput;
  }
  GaussianTrackDensity::State densityState(trackVector.size());

  // Calculate z seed position
  auto zAndWidthRes = m_cfg.trackDensityEstimator.globalMaximumWithWidth(
      densityState, trackVector);
  if (!zAndWidthRes.ok()) {
    return zAndWidthRes.error();
  }
  const auto& zAndWidthOpt = *zAndWidthRes;
  if (!zAndWidthOpt) {
    return std::vector<Vertex>();
  }
  const auto& zAndWidth = *zAndWidthOpt;

  double z = zAndWidth.first;

  // Calculate seed position
  // The default reference has zero position, time and covariance.
  const Vertex reference = vertexingOptions.seedReference.value_or(Vertex{});
  Vector4 seedPos = reference.fullPosition() + Vector4(0., 0., z, 0.);

  Vertex returnVertex = Vertex(seedPos);

  SquareMatrix4 seedCov = reference.fullCovariance();

  // Update the reference's longitudinal uncertainty from the peak width.
  if (seedCov != SquareMatrix4::Zero() && std::isnormal(zAndWidth.second)) {
    seedCov(eZ, eZ) = zAndWidth.second * zAndWidth.second;
  }

  returnVertex.setFullCovariance(seedCov);

  return std::vector<Vertex>{returnVertex};
}
