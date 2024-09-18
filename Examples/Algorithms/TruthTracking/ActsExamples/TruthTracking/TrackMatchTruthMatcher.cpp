// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackMatchTruthMatcher.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace ActsExamples {

TrackMatchTruthMatcher::TrackMatchTruthMatcher(const Config& config,
                                     Acts::Logging::Level level)
    : IAlgorithm("TrackMatchTruthMatcher", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles");
  }
  if (m_cfg.inputMeasurementParticlesMapVT.empty()) {
    throw std::invalid_argument("Missing input measurement particles map");
  }
  if (m_cfg.inputMeasurementParticlesMapMS.empty()) {
    throw std::invalid_argument("Missing input measurement particles map");
  }

  m_inputTracksVT.initialize(m_cfg.inputTracksVT);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurementParticlesMapVT.initialize(m_cfg.inputMeasurementParticlesMapVT);
  m_inputTracksMS.initialize(m_cfg.inputTracksMS);
  m_inputMeasurementParticlesMapMS.initialize(m_cfg.inputMeasurementParticlesMapMS);
}

ActsExamples::ProcessCode TrackMatchTruthMatcher::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input tracks
  const auto& tracksMS = m_inputTracksMS(ctx);
  const auto& tracksVT = m_inputTracksVT(ctx);

  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& hitParticlesMapMS = m_inputMeasurementParticlesMapMS(ctx);
  const auto& hitParticlesMapVT = m_inputMeasurementParticlesMapVT(ctx);

  // TODO this may be computed in a separate algorithm
  // TODO can we wire this through?
  std::map<SimBarcode, std::size_t> particleTruthHitCount;
  for (const auto& [_, pid] : hitParticlesMap) {
    particleTruthHitCount[pid]++;
  }

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  for (const auto& track : tracksVT) {
    // Get the majority truth particle to this track
    identifyContributingParticles(hitParticlesMap, track, particleHitCounts);
    if (particleHitCounts.empty()) {
      ACTS_DEBUG(
          "No truth particle associated with this trajectory with tip index = "
          << track.tipIndex());
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
        static_cast<double>(nMajorityHits) / track.nMeasurements() >=
        m_cfg.matchingRatio;
    const bool truthMatched =
        static_cast<double>(nMajorityHits) /
            particleTruthHitCount.at(majorityParticleId) >=
        m_cfg.matchingRatio;

    if ((!m_cfg.doubleMatching && recoMatched) ||
        (m_cfg.doubleMatching && recoMatched && truthMatched)) {
      auto& trackParticleMatch = trackParticleMatching[track.index()] = {
          TrackMatchClassification::Matched, majorityParticleId,
          particleHitCounts};

      auto& particleTrackMatch = particleTrackMatching[majorityParticleId];
      if (!particleTrackMatch.track) {
        particleTrackMatch.track = track.index();
      } else {
        // we already have a track associated with this particle and have to
        // resolve the ambiguity.
        // we will use the track with more hits and smaller chi2
        const auto& otherTrack =
            tracks.getTrack(particleTrackMatch.track.value());
        if (otherTrack.nMeasurements() < track.nMeasurements() ||
            otherTrack.chi2() > track.chi2()) {
          trackParticleMatching[otherTrack.index()].classification =
              TrackMatchClassification::Duplicate;
          particleTrackMatch.track = track.index();
        } else {
          trackParticleMatch.classification =
              TrackMatchClassification::Duplicate;
        }

        ++particleTrackMatch.duplicates;
      }
    } else {
      trackParticleMatching[track.index()] = {TrackMatchClassification::Fake,
                                              std::nullopt, particleHitCounts};

      auto& particleTrackMatch = particleTrackMatching[majorityParticleId];
      ++particleTrackMatch.fakes;
    }
  }


  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
