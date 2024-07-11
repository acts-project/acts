// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TruthVertexSeeder.hpp"

#include "VertexingHelpers.hpp"

namespace ActsExamples {

TruthVertexSeeder::TruthVertexSeeder(const Config &cfg) : m_cfg(cfg) {}

Acts::Result<std::vector<Acts::Vertex>> TruthVertexSeeder::find(
    const std::vector<Acts::InputTrack> & /*trackVector*/,
    const Acts::VertexingOptions &vertexingOptions,
    Acts::IVertexFinder::State &anyState) const {
  auto &state = anyState.template as<State>();

  if (state.nextVertexIndex >= m_cfg.vertices.size()) {
    std::cout << "TruthVertexSeeder: no more vertices to seed" << std::endl;
    return std::vector<Acts::Vertex>();
  }

  const auto &truthVertex = m_cfg.vertices[state.nextVertexIndex];
  ++state.nextVertexIndex;

  Acts::Vertex converted;
  converted.fullPosition().z() = truthVertex.position().z();
  if (m_cfg.useXY) {
    converted.fullPosition().x() = truthVertex.position().x();
    converted.fullPosition().y() = truthVertex.position().y();
  }
  if (m_cfg.useTime) {
    converted.setTime(truthVertex.time());
  }

  Acts::SquareMatrix4 seedCov = vertexingOptions.constraint.fullCovariance();
  converted.setFullCovariance(seedCov);

  std::vector<Acts::Vertex> vertices;
  vertices.push_back(converted);
  std::cout << "TruthVertexSeeder: seeded vertex at "
            << converted.position().transpose() << std::endl;
  return vertices;
}

Acts::IVertexFinder::State TruthVertexSeeder::makeState(
    const Acts::MagneticFieldContext & /*mctx*/) const {
  return Acts::IVertexFinder::State{State{}};
}

void TruthVertexSeeder::setTracksToRemove(
    Acts::IVertexFinder::State & /*anyState*/,
    const std::vector<Acts::InputTrack> & /*removedTracks*/) const {
  // nothing to do
}

}  // namespace ActsExamples
