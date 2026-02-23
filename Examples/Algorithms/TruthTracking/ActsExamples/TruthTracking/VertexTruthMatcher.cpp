// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/VertexTruthMatcher.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Utilities/VertexTruthUtility.hpp"

#include <optional>
#include <stdexcept>

namespace ActsExamples {

VertexTruthMatcher::VertexTruthMatcher(const Config& config,
                                       Acts::Logging::Level level)
    : IAlgorithm("VertexTruthMatcher",
                 Acts::getDefaultLogger("VertexTruthMatcher", level)),
      m_cfg(config) {
  if (m_cfg.inputVertices.empty()) {
    throw std::invalid_argument("Collection with vertices missing");
  }
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Collection with tracks missing");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Collection with particles missing");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument(
        "Collection with track-particle matching missing");
  }
  if (m_cfg.outputVertexTruthMatching.empty()) {
    throw std::invalid_argument("Missing output vertex-truth matching");
  }
  if (m_cfg.outputTruthVertexMatching.empty()) {
    throw std::invalid_argument("Missing output truth-vertex matching");
  }

  m_inputVertices.initialize(m_cfg.inputVertices);
  m_inputTracks.maybeInitialize(m_cfg.inputTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.maybeInitialize(
      m_cfg.inputTrackParticleMatching);
  m_outputVertexTruthMatching.initialize(m_cfg.outputVertexTruthMatching);
  m_outputTruthVertexMatching.initialize(m_cfg.outputTruthVertexMatching);
}

ActsExamples::ProcessCode VertexTruthMatcher::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input reconstructed vertices
  const auto& vertices = m_inputVertices(ctx);
  // Read input tracks
  const auto& tracks = m_inputTracks(ctx);
  // Read truth particle input collection
  const SimParticleContainer& particles = m_inputParticles(ctx);
  // Read track-particle matching
  const TrackParticleMatching& trackParticleMatching =
      m_inputTrackParticleMatching(ctx);

  std::vector<VertexToTruthMatching> recoToTruthMatching;
  std::map<SimVertexBarcode, VertexToRecoMatching> truthToRecoMatching;

  // Do truth matching for each reconstructed vertex
  for (const auto& [vtxIndex, vtx] : Acts::enumerate(vertices)) {
    // Reconstructed tracks that contribute to the reconstructed vertex
    const std::vector<Acts::TrackAtVertex>& tracksAtVtx = vtx.tracks();

    // Containers for storing truth particles and truth vertices that
    // contribute to the reconstructed vertex
    std::vector<std::pair<SimVertexBarcode, double>> contributingTruthVertices;

    double totalTrackWeight = 0;
    for (const Acts::TrackAtVertex& trk : tracksAtVtx) {
      if (trk.trackWeight < m_cfg.minTrkWeight) {
        continue;
      }

      totalTrackWeight += trk.trackWeight;

      std::optional<ConstTrackProxy> trackOpt = findTrack(tracks, trk);
      if (!trackOpt.has_value()) {
        ACTS_DEBUG("Track has no matching input track.");
        continue;
      }
      const ConstTrackProxy& inputTrk = *trackOpt;
      const SimParticle* particle =
          findParticle(particles, trackParticleMatching, inputTrk, logger());
      if (particle == nullptr) {
        ACTS_VERBOSE("Track has no matching truth particle.");
      } else {
        contributingTruthVertices.emplace_back(
            SimBarcode{particle->particleId()}.vertexId(), trk.trackWeight);
      }
    }

    // Find true vertex that contributes most to the reconstructed vertex
    std::map<SimVertexBarcode, std::pair<int, double>> fmap;
    for (const auto& [vtxId, weight] : contributingTruthVertices) {
      ++fmap[vtxId].first;
      fmap[vtxId].second += weight;
    }
    double truthMajorityVertexTrackWeights = 0;
    std::optional<SimVertexBarcode> truthMajorityVertexId = std::nullopt;
    for (const auto& [vtxId, counter] : fmap) {
      if (counter.second > truthMajorityVertexTrackWeights) {
        truthMajorityVertexId = vtxId;
        truthMajorityVertexTrackWeights = counter.second;
      }
    }

    double sumPt2 = calcSumPt2(vtx, m_cfg.minTrkWeight);

    double vertexMatchFraction =
        (totalTrackWeight > 0
             ? truthMajorityVertexTrackWeights / totalTrackWeight
             : 0.);
    RecoVertexClassification recoVertexClassification =
        RecoVertexClassification::Unknown;

    if (vertexMatchFraction >= m_cfg.vertexMatchThreshold) {
      recoVertexClassification = RecoVertexClassification::Clean;
    } else {
      recoVertexClassification = RecoVertexClassification::Merged;
    }

    recoToTruthMatching.push_back({truthMajorityVertexId, totalTrackWeight,
                                   truthMajorityVertexTrackWeights,
                                   vertexMatchFraction,
                                   recoVertexClassification});

    auto& recoToTruth = recoToTruthMatching.back();

    // We have to decide if this reco vertex is a split vertex.
    if (!truthMajorityVertexId.has_value()) {
      // No truth vertex matched to this reconstructed vertex
      ACTS_DEBUG("No truth vertex matched to this reconstructed vertex.");
      continue;
    }

    if (auto it = truthToRecoMatching.find(truthMajorityVertexId.value());
        it != truthToRecoMatching.end()) {
      // This truth vertex is already matched to a reconstructed vertex so we
      // are dealing with a split vertex.

      // We have to decide which of the two reconstructed vertices is the
      // split vertex.
      if (sumPt2 <= it->second.recoSumPt2) {
        // Since the sumPt2 is smaller we can simply call this a split vertex

        recoToTruth.classification = RecoVertexClassification::Split;

        // Keep the existing truth to reco matching
      } else {
        // The sumPt2 is larger, so we call the other vertex a split vertex.

        auto& otherRecoToTruth = recoToTruthMatching.at(it->second.recoIndex);
        // Swap the classification
        recoToTruth.classification = otherRecoToTruth.classification;
        otherRecoToTruth.classification = RecoVertexClassification::Split;

        // Overwrite the truth to reco matching
        it->second = {vtxIndex, sumPt2};
      }
    } else {
      truthToRecoMatching[truthMajorityVertexId.value()] = {vtxIndex, sumPt2};
    }
  }

  m_outputVertexTruthMatching(ctx, std::move(recoToTruthMatching));
  m_outputTruthVertexMatching(ctx, std::move(truthToRecoMatching));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
