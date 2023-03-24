// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthVertexFinder.hpp"

#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

ActsExamples::TruthVertexFinder::TruthVertexFinder(const Config& config,
                                                   Acts::Logging::Level level)
    : IAlgorithm("TruthVertexFinder", level), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputProtoVertices.initialize(m_cfg.outputProtoVertices);
}

ActsExamples::ProcessCode ActsExamples::TruthVertexFinder::execute(
    const AlgorithmContext& ctx) const {
  // prepare input and output collections
  ACTS_VERBOSE("Reading particles from " << m_cfg.inputParticles);
  const auto& particles = m_inputParticles(ctx);
  ProtoVertexContainer protoVertices;
  ACTS_VERBOSE("Have " << particles.size() << " particles");

  // assumes the begin/end iterator references the particles container
  auto addProtoVertex = [&](SimParticleContainer::const_iterator begin,
                            const SimParticleContainer::const_iterator& end) {
    ProtoVertex protoVertex;
    protoVertex.reserve(std::distance(begin, end));
    // determine each particle index
    for (; begin != end; ++begin) {
      protoVertex.push_back(std::distance(particles.begin(), begin));
    }
    protoVertices.push_back(std::move(protoVertex));
  };

  if (m_cfg.excludeSecondaries) {
    // if secondaries are excluded, the `separateSecondaries` flag has no effect
    // since there will be no secondary vertices to separate
    for (auto&& [vtxId, vtxParticles] : groupBySecondaryVertex(particles)) {
      if (vtxId.vertexSecondary() != 0u) {
        continue;
      }
      addProtoVertex(vtxParticles.begin(), vtxParticles.end());
    }
  } else {
    // particles from secondary vertices should be included
    if (m_cfg.separateSecondaries) {
      // secondary particles are added to separate secondary vertices
      for (auto&& [vtxId, vtxParticles] : groupBySecondaryVertex(particles)) {
        addProtoVertex(vtxParticles.begin(), vtxParticles.end());
      }
    } else {
      // secondary particles are included in the primary vertex
      for (auto&& [vtxId, vtxParticles] : groupByPrimaryVertex(particles)) {
        addProtoVertex(vtxParticles.begin(), vtxParticles.end());
      }
    }
  }

  ACTS_VERBOSE("Write " << protoVertices.size() << " proto vertex to "
                        << m_cfg.outputProtoVertices);
  m_outputProtoVertices(ctx, std::move(protoVertices));
  return ProcessCode::SUCCESS;
}
