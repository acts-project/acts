// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/ProtoTrackTruthMatcher.hpp"

#include "Acts/Utilities/Enumerate.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace ActsExamples {

ProtoTrackTruthMatcher::ProtoTrackTruthMatcher(const Config& config,
                                               Acts::Logging::Level level)
    : IAlgorithm("ProtoTrackTruthMatcher", level), m_cfg(config) {
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing input proto tracks");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input measurement particles map");
  }
  if (m_cfg.outputProtoTrackParticleMatching.empty()) {
    throw std::invalid_argument(
        "Missing output proto track particles matching");
  }
  if (m_cfg.outputParticleProtoTrackMatching.empty()) {
    throw std::invalid_argument("Missing output particle proto track matching");
  }

  m_inputProtoTracks.initialize(m_cfg.inputProtoTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_outputProtoTrackParticleMatching.initialize(
      m_cfg.outputProtoTrackParticleMatching);
  m_outputParticleProtoTrackMatching.initialize(
      m_cfg.outputParticleProtoTrackMatching);
}

ActsExamples::ProcessCode ProtoTrackTruthMatcher::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input proto tracks
  const auto& protoTracks = m_inputProtoTracks(ctx);

  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);

  ProtoTrackParticleMatching protoTrackParticleMatching;
  ParticleProtoTrackMatching particleProtoTrackMatching;

  // TODO this may be computed in a separate algorithm
  // TODO can we wire this through?
  std::map<SimBarcode, std::size_t> particleTruthHitCount;
  for (const auto& [_, pid] : hitParticlesMap) {
    particleTruthHitCount[pid]++;
  }

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  for (const auto& [index, protoTrack] : Acts::enumerate(protoTracks)) {
    // Get the majority truth particle to this track
    identifyContributingParticles(hitParticlesMap, protoTrack,
                                  particleHitCounts);
    if (particleHitCounts.empty()) {
      ACTS_DEBUG("No truth particle associated with this proto track "
                 << index);
      continue;
    }

    // Get the majority particleId and majority particle counts
    // Note that the majority particle might not be in the truth seeds
    // collection
    ActsFatras::Barcode majorityParticleId =
        particleHitCounts.front().particleId;
    std::size_t nMajorityHits = particleHitCounts.front().hitCount;

    if (particles.find(majorityParticleId) == particles.end()) {
      ACTS_DEBUG(
          "The majority particle is not in the input particle collection, "
          "majorityParticleId = "
          << majorityParticleId);
      continue;
    }

    // Check if the trajectory is matched with truth.
    // If not, it will be classified as 'fake'
    const bool recoMatched =
        static_cast<double>(nMajorityHits) / protoTrack.size() >=
        m_cfg.matchingRatio;
    const bool truthMatched =
        static_cast<double>(nMajorityHits) /
            particleTruthHitCount.at(majorityParticleId) >=
        m_cfg.matchingRatio;

    if ((!m_cfg.doubleMatching && recoMatched) ||
        (m_cfg.doubleMatching && recoMatched && truthMatched)) {
      auto& trackParticleMatch = protoTrackParticleMatching[index] = {
          TrackMatchClassification::Matched, majorityParticleId,
          particleHitCounts};

      auto& particleTrackMatch = particleProtoTrackMatching[majorityParticleId];
      if (!particleTrackMatch.track) {
        particleTrackMatch.track = index;
      } else {
        // we already have a track associated with this particle and have to
        // resolve the ambiguity.
        // we will use the track with more hits and smaller chi2
        const auto& otherProtoTrack =
            protoTracks.at(particleTrackMatch.track.value());
        if (otherProtoTrack.size() < protoTrack.size()) {
          protoTrackParticleMatching[particleTrackMatch.track.value()]
              .classification = TrackMatchClassification::Duplicate;
          particleTrackMatch.track = index;
        } else {
          trackParticleMatch.classification =
              TrackMatchClassification::Duplicate;
        }

        ++particleTrackMatch.duplicates;
      }
    } else {
      protoTrackParticleMatching[index] = {TrackMatchClassification::Fake,
                                           std::nullopt, particleHitCounts};

      auto& particleTrackMatch = particleProtoTrackMatching[majorityParticleId];
      ++particleTrackMatch.fakes;
    }
  }

  m_outputProtoTrackParticleMatching(ctx,
                                     std::move(protoTrackParticleMatching));
  m_outputParticleProtoTrackMatching(ctx,
                                     std::move(particleProtoTrackMatching));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
