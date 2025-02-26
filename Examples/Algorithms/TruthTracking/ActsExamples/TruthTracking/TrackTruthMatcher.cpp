// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackTruthMatcher.hpp"

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace ActsExamples {

TrackTruthMatcher::TrackTruthMatcher(const Config& config,
                                     Acts::Logging::Level level)
    : IAlgorithm("TrackTruthMatcher", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input measurement particles map");
  }
  if (m_cfg.outputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing output track particles matching");
  }
  if (m_cfg.outputParticleTrackMatching.empty()) {
    throw std::invalid_argument("Missing output particle track matching");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_outputTrackParticleMatching.initialize(m_cfg.outputTrackParticleMatching);
  m_outputParticleTrackMatching.initialize(m_cfg.outputParticleTrackMatching);
}

ActsExamples::ProcessCode TrackTruthMatcher::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input tracks
  const auto& tracks = m_inputTracks(ctx);

  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);

  TrackParticleMatching trackParticleMatching;
  ParticleTrackMatching particleTrackMatching;

  // TODO this may be computed in a separate algorithm
  // TODO can we wire this through?
  std::map<SimBarcode, std::size_t> particleTruthHitCount;
  for (const auto& [_, pid] : hitParticlesMap) {
    particleTruthHitCount[pid]++;
  }

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;

  std::size_t failed = 0;
  for (const auto& track : tracks) {
    // Get the majority truth particle to this track
    identifyContributingParticles(hitParticlesMap, track, particleHitCounts);
    if (particleHitCounts.empty()) {
      ACTS_VERBOSE(
          "No truth particle associated with this trajectory with tip index = "
          << track.tipIndex());
      failed++;
      continue;
    }

    // Get the majority particleId and majority particle counts
    // Note that the majority particle might not be in the truth seeds
    // collection
    ActsFatras::Barcode majorityParticleId =
        particleHitCounts.front().particleId;
    std::size_t nMajorityHits = particleHitCounts.front().hitCount;

    if (!particles.contains(majorityParticleId)) {
      ACTS_VERBOSE(
          "The majority particle is not in the input particle collection, "
          "majorityParticleId = "
          << majorityParticleId);
      failed++;
      continue;
    }

    // For now only do the weighted for reco matching
    auto weightedRecoMatching = [&]() {
      double denominator = 0.0;
      std::map<ActsFatras::Barcode, double> nominators;
      std::size_t nReweighted = 0;
      for (const auto& state : track.trackStatesReversed()) {
        // no truth info with non-measurement state
        if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          continue;
        }
        // register all particles that generated this hit
        IndexSourceLink sl =
            state.getUncalibratedSourceLink().template get<IndexSourceLink>();
        auto vol = sl.geometryId().volume();

        double weight = m_cfg.reweightVolumes.contains(vol)
                            ? m_cfg.reweightVolumes.at(vol)
                            : 1.0;
        denominator += weight;

        nReweighted += static_cast<int>(weight == 1.0);

        auto hitIndex = sl.index();
        for (auto [_, particleId] :
             makeRange(hitParticlesMap.equal_range(hitIndex))) {
          nominators[particleId] += weight;
        }
      }
      auto nominator =
          std::ranges::max_element(nominators, std::less{}, [](const auto& a) {
            return a.second;
          })->second;
      return nominator / denominator;
    };

    // Check if the trajectory is matched with truth.
    // If not, it will be classified as 'fake'
    const double recoMatchingProb =
        m_cfg.reweightVolumes.empty()
            ? static_cast<double>(nMajorityHits) / track.nMeasurements()
            : weightedRecoMatching();

    /*auto p1 = weightedRecoMatching();
    auto p2 = static_cast<double>(nMajorityHits) / track.nMeasurements();
    if( !m_cfg.reweightVolumes.empty() && p1 != p2) {
      ACTS_INFO("reweight prob: " << p1 << ", unweighted: " << p2);
    }*/

    const bool recoMatched = recoMatchingProb >= m_cfg.matchingRatio;

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

    trackParticleMatching.at(track.index()).majorityMatchingProbability =
        recoMatchingProb;
  }

  ACTS_DEBUG("Matching failed for " << failed << " tracks");

  m_outputTrackParticleMatching(ctx, std::move(trackParticleMatching));
  m_outputParticleTrackMatching(ctx, std::move(particleTrackMatching));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
