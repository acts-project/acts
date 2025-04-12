// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthVertexFinder.hpp"

#include "ActsExamples/EventData/ProtoVertex.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ActsExamples {

namespace {

using ConstTrackProxy = ConstTrackContainer::ConstTrackProxy;
using TrackIndex = ConstTrackProxy::IndexType;

struct TrackMatchEntry {
  std::optional<SimBarcode> particle;

  /// Number of hits on the track that are associated to a particle
  /// Sorted by decreasing number of hits
  std::vector<ParticleHitCount> contributingParticles;
};

struct ParticleMatchEntry {
  std::optional<TrackIndex> track;
  std::uint32_t duplicates{};
  std::uint32_t fakes{};
};

}  // namespace

TruthVertexFinder::TruthVertexFinder(const Config& config,
                                     Acts::Logging::Level level)
    : IAlgorithm("TruthVertexFinder", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputParticleTrackMatching.empty()) {
    throw std::invalid_argument("Missing input particle-track matching");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputParticleTrackMatching.initialize(m_cfg.inputParticleTrackMatching);
  m_outputProtoVertices.initialize(m_cfg.outputProtoVertices);
}

ProcessCode TruthVertexFinder::execute(const AlgorithmContext& ctx) const {
  // prepare input and output collections
  const auto& tracks = m_inputTracks(ctx);
  const auto& particles = m_inputParticles(ctx);
  const auto& particleTrackMatching = m_inputParticleTrackMatching(ctx);

  ACTS_VERBOSE("Have " << particles.size() << " particles");
  ACTS_VERBOSE("Have " << tracks.size() << " tracks");

  std::unordered_map<SimBarcode, std::vector<TrackIndex>> protoVertexTrackMap;

  for (const auto& [particle, trackMatch] : particleTrackMatching) {
    auto vtxId = SimBarcode(particle).setParticle(0).setSubParticle(0);

    if (trackMatch.track) {
      protoVertexTrackMap[vtxId].push_back(trackMatch.track.value());
    }
  }

  ProtoVertexContainer protoVertices;

  // assumes the begin/end iterator references the particles container
  auto addProtoVertex = [&](const std::vector<TrackIndex>& vtxTracks) {
    ProtoVertex protoVertex;
    protoVertex.reserve(vtxTracks.size());
    for (const auto& track : vtxTracks) {
      protoVertex.push_back(track);
    }
    protoVertices.push_back(std::move(protoVertex));
  };

  if (m_cfg.excludeSecondaries) {
    // if secondaries are excluded, the `separateSecondaries` flag has no effect
    // since there will be no secondary vertices to separate
    for (auto&& [vtxId, vtxTracks] : protoVertexTrackMap) {
      if (vtxId.vertexSecondary() != 0u) {
        continue;
      }
      addProtoVertex(vtxTracks);
    }
  } else {
    // particles from secondary vertices should be included
    if (m_cfg.separateSecondaries) {
      // secondary particles are added to separate secondary vertices
      for (auto&& [vtxId, vtxTracks] : protoVertexTrackMap) {
        addProtoVertex(vtxTracks);
      }
    } else {
      // secondary particles are included in the primary vertex

      std::unordered_map<SimBarcode, std::vector<TrackIndex>>
          protoVertexTrackMap2;
      for (auto&& [vtxId, vtxTracks] : protoVertexTrackMap) {
        auto vtxId2 = SimBarcode(vtxId).setVertexSecondary(0);
        protoVertexTrackMap2[vtxId2].insert(protoVertexTrackMap2[vtxId].end(),
                                            vtxTracks.begin(), vtxTracks.end());
      }

      for (auto&& [vtxId, vtxTracks] : protoVertexTrackMap2) {
        addProtoVertex(vtxTracks);
      }
    }
  }

  ACTS_VERBOSE("Write " << protoVertices.size() << " proto vertex to "
                        << m_cfg.outputProtoVertices);

  m_outputProtoVertices(ctx, std::move(protoVertices));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
