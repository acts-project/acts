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
#include "ActsExamples/EventData/SimVertex.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
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
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing input track-particle matching");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
  m_outputProtoVertices.initialize(m_cfg.outputProtoVertices);
}

ProcessCode TruthVertexFinder::execute(const AlgorithmContext& ctx) const {
  // prepare input and output collections
  const auto& tracks = m_inputTracks(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);

  ACTS_DEBUG("Have " << tracks.size() << " tracks");

  // first step is to sort the tracks into potential proto vertices using the
  // track-particle matching information.
  // at this stage we use the full vertex ID, including the secondary vertex ID.

  std::unordered_map<SimVertexBarcode, std::vector<TrackIndex>>
      protoVertexTrackMap;

  for (const auto& track : tracks) {
    auto trackMatchIt = trackParticleMatching.find(track.index());
    if (trackMatchIt == trackParticleMatching.end()) {
      continue;
    }
    const auto& trackMatch = trackMatchIt->second;

    // get the particle associated to the track
    auto particleOpt = trackMatchIt->second.particle;
    if (!particleOpt) {
      continue;
    }
    auto barcode = *particleOpt;

    // Skip fake and duplicate tracks
    if (trackMatch.classification != TrackMatchClassification::Matched) {
      continue;
    }

    // derive the vertex ID from the barcode
    SimVertexBarcode vertexId = SimVertexBarcode(barcode);

    // add the track to the proto vertex map
    protoVertexTrackMap[vertexId].push_back(track.index());
  }

  // in the second step we separate secondary vertices based on the
  // configuration.

  ProtoVertexContainer protoVertices;

  // assumes the begin/end iterator references the particles container
  auto addProtoVertex = [&](const std::vector<TrackIndex>& vertexTracks) {
    ProtoVertex protoVertex;
    protoVertex.reserve(vertexTracks.size());
    for (const auto& track : vertexTracks) {
      protoVertex.push_back(track);
    }
    protoVertices.push_back(std::move(protoVertex));
  };

  if (m_cfg.excludeSecondaries) {
    // if secondaries are excluded, the `separateSecondaries` flag has no effect
    // since there will be no secondary vertices to separate
    for (auto&& [vertexId, vertexTracks] : protoVertexTrackMap) {
      if (vertexId.vertexSecondary() != 0u) {
        continue;
      }
      addProtoVertex(vertexTracks);
    }
  } else {
    // particles from secondary vertices should be included
    if (m_cfg.separateSecondaries) {
      // secondary particles are added to separate secondary vertices
      for (auto&& [vertexId, vertexTracks] : protoVertexTrackMap) {
        addProtoVertex(vertexTracks);
      }
    } else {
      // secondary particles are included in the primary vertex

      std::unordered_map<SimVertexBarcode, std::vector<TrackIndex>>
          protoVertexTrackMap2;
      for (auto&& [vertexId, vertexTracks] : protoVertexTrackMap) {
        auto vertexId2 = SimVertexBarcode(vertexId).setVertexSecondary(0);
        auto& vertexTracks2 = protoVertexTrackMap2[vertexId2];
        vertexTracks2.insert(vertexTracks2.end(), vertexTracks.begin(),
                             vertexTracks.end());
      }

      for (auto&& [vertexId, vertexTracks] : protoVertexTrackMap2) {
        addProtoVertex(vertexTracks);
      }
    }
  }

  ACTS_DEBUG("Write " << protoVertices.size() << " proto vertex to "
                      << m_cfg.outputProtoVertices);

  m_outputProtoVertices(ctx, std::move(protoVertices));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
