// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "TruthVertexSeeder.hpp"

namespace ActsExamples {

TruthVertexSeeder::TruthVertexSeeder(const Config &cfg) : m_cfg(cfg) {}

Acts::Result<VertexContainer> TruthVertexSeeder::find(
    const std::vector<Acts::InputTrack> & /*trackVector*/,
    const Acts::VertexingOptions &vertexingOptions,
    Acts::IVertexFinder::State &anyState) const {
  auto &state = anyState.template as<State>();

  if (state.nextVertexIndex >= state.truthVertices.size()) {
    return VertexContainer();
  }

  VertexContainer seeds;
  for (std::size_t i = 0; i < m_cfg.simultaneousSeeds; ++i) {
    if (state.nextVertexIndex >= state.truthVertices.size()) {
      break;
    }

    const auto &truthVertex = state.truthVertices[state.nextVertexIndex];
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

    seeds.push_back(converted);
  }

  return seeds;
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
