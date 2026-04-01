// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackFinderPerformanceCollector.hpp"

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <utility>

namespace ActsExamples {

TrackFinderPerformanceCollector::TrackFinderPerformanceCollector(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(std::move(cfg)),
      m_logger(std::move(logger)),
      m_effPlotTool(m_cfg.effPlotToolConfig, m_logger->level()),
      m_fakePlotTool(m_cfg.fakePlotToolConfig, m_logger->level()),
      m_duplicationPlotTool(m_cfg.duplicationPlotToolConfig, m_logger->level()),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig,
                             m_logger->level()),
      m_trackQualityPlotTool(m_cfg.trackQualityPlotToolConfig,
                             m_logger->level()) {
  for (const auto& [key, _] : m_cfg.subDetectorTrackSummaryVolumes) {
    TrackSummaryPlotTool::Config subConfig = m_cfg.trackSummaryPlotToolConfig;
    subConfig.prefix = key;
    m_subDetectorSummaryTools.emplace(
        std::piecewise_construct, std::forward_as_tuple(key),
        std::forward_as_tuple(subConfig, m_logger->level()));
  }
}

ProcessCode TrackFinderPerformanceCollector::fill(
    const AlgorithmContext& ctx, const ConstTrackContainer& tracks,
    const SimParticleContainer& particles,
    const TrackParticleMatching& trackParticleMatching,
    const ParticleTrackMatching& particleTrackMatching,
    const InverseMultimap<SimBarcode>& particleMeasurementsMap) {
  std::size_t unmatched = 0;
  std::size_t missingRefSurface = 0;

  for (const auto& track : tracks) {
    m_stats.nTotalTracks++;

    if (!track.hasReferenceSurface()) {
      missingRefSurface++;
      continue;
    }

    Acts::BoundTrackParameters fittedParameters =
        track.createParametersAtReference();

    m_trackSummaryPlotTool.fill(fittedParameters, track.nTrackStates(),
                                track.nMeasurements(), track.nOutliers(),
                                track.nHoles(), track.nSharedHits());

    for (const auto& [key, volumes] : m_cfg.subDetectorTrackSummaryVolumes) {
      std::size_t nTrackStates{};
      std::size_t nMeasurements{};
      std::size_t nOutliers{};
      std::size_t nHoles{};
      std::size_t nSharedHits{};

      for (auto state : track.trackStatesReversed()) {
        if (!state.hasReferenceSurface() ||
            !volumes.contains(state.referenceSurface().geometryId().volume())) {
          continue;
        }
        nTrackStates++;
        nMeasurements +=
            static_cast<std::size_t>(state.typeFlags().isMeasurement());
        nOutliers += static_cast<std::size_t>(state.typeFlags().isOutlier());
        nHoles += static_cast<std::size_t>(state.typeFlags().isHole());
        nSharedHits +=
            static_cast<std::size_t>(state.typeFlags().isSharedHit());
      }
      m_subDetectorSummaryTools.at(key).fill(fittedParameters, nTrackStates,
                                             nMeasurements, nOutliers, nHoles,
                                             nSharedHits);
    }

    auto imatched = trackParticleMatching.find(track.index());
    if (imatched == trackParticleMatching.end()) {
      unmatched++;
      continue;
    }

    const auto& particleMatch = imatched->second;

    if (particleMatch.classification == TrackMatchClassification::Fake) {
      m_stats.nTotalFakeTracks++;
    }
    if (particleMatch.classification == TrackMatchClassification::Duplicate) {
      m_stats.nTotalDuplicateTracks++;
    }

    m_fakePlotTool.fill(fittedParameters, particleMatch.classification ==
                                              TrackMatchClassification::Fake);
    m_duplicationPlotTool.fill(
        fittedParameters,
        particleMatch.classification == TrackMatchClassification::Duplicate);

    if (particleMatch.particle.has_value() &&
        particleMeasurementsMap.contains(particleMatch.particle.value())) {
      const auto measurements =
          particleMeasurementsMap.equal_range(particleMatch.particle.value());

      std::size_t nTrackMeasurements =
          track.nMeasurements() + track.nOutliers();
      std::size_t nMatchedHits =
          particleMatch.contributingParticles.front().hitCount;
      std::size_t nParticleHits =
          std::distance(measurements.first, measurements.second);

      double completeness = static_cast<double>(nMatchedHits) / nParticleHits;
      double purity = static_cast<double>(nMatchedHits) / nTrackMeasurements;

      m_trackQualityPlotTool.fill(fittedParameters, completeness, purity);
    }
  }

  if (unmatched > 0) {
    ACTS_VERBOSE("No matching information found for " << unmatched
                                                      << " tracks");
  }
  if (missingRefSurface > 0) {
    ACTS_VERBOSE("Reference surface was missing for " << missingRefSurface
                                                      << " tracks");
  }

  for (const auto& particle : particles) {
    auto particleId = particle.particleId();

    std::size_t nMatchedTracks = 0;
    std::size_t nFakeTracks = 0;
    bool isReconstructed = false;
    if (auto imatched = particleTrackMatching.find(particleId);
        imatched != particleTrackMatching.end()) {
      isReconstructed = imatched->second.track.has_value();
      nMatchedTracks = (isReconstructed ? 1 : 0) + imatched->second.duplicates;

      m_stats.nTotalMatchedTracks += nMatchedTracks;
      m_stats.nTotalMatchedParticles += isReconstructed ? 1 : 0;

      if (nMatchedTracks > 1) {
        m_stats.nTotalDuplicateParticles += 1;
      }

      nFakeTracks = imatched->second.fakes;
      if (nFakeTracks > 0) {
        m_stats.nTotalFakeParticles += 1;
      }
    }

    double minDeltaR = -1;
    for (const auto& closeParticle : particles) {
      if (closeParticle.particleId() == particleId) {
        continue;
      }
      double distance = Acts::VectorHelpers::deltaR(particle.direction(),
                                                    closeParticle.direction());
      if (minDeltaR == -1 || distance < minDeltaR) {
        minDeltaR = distance;
      }
    }

    m_effPlotTool.fill(ctx.geoContext, particle.initialState(), minDeltaR,
                       isReconstructed);
    m_duplicationPlotTool.fill(particle.initialState(), nMatchedTracks);
    m_fakePlotTool.fill(particle.initialState(), nMatchedTracks, nFakeTracks);

    m_stats.nTotalParticles += 1;
  }

  return ProcessCode::SUCCESS;
}

void TrackFinderPerformanceCollector::logSummary() const {
  const Acts::Logger& log = *m_logger;
  float eff_tracks =
      static_cast<float>(m_stats.nTotalMatchedTracks) / m_stats.nTotalTracks;
  float fakeRatio_tracks =
      static_cast<float>(m_stats.nTotalFakeTracks) / m_stats.nTotalTracks;
  float duplicationRatio_tracks =
      static_cast<float>(m_stats.nTotalDuplicateTracks) / m_stats.nTotalTracks;

  float eff_particle = static_cast<float>(m_stats.nTotalMatchedParticles) /
                       m_stats.nTotalParticles;
  float fakeRatio_particle =
      static_cast<float>(m_stats.nTotalFakeParticles) / m_stats.nTotalParticles;
  float duplicationRatio_particle =
      static_cast<float>(m_stats.nTotalDuplicateParticles) /
      m_stats.nTotalParticles;

  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::DEBUG,
      "nTotalTracks                = " << m_stats.nTotalTracks);
  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::DEBUG,
      "nTotalMatchedTracks         = " << m_stats.nTotalMatchedTracks);
  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::DEBUG,
      "nTotalDuplicateTracks       = " << m_stats.nTotalDuplicateTracks);
  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::DEBUG,
      "nTotalFakeTracks            = " << m_stats.nTotalFakeTracks);

  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::INFO,
      "Efficiency with tracks (nMatchedTracks/ nAllTracks) = " << eff_tracks);
  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::INFO,
      "Fake ratio with tracks (nFakeTracks/nAllTracks) = " << fakeRatio_tracks);
  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::INFO,
      "Duplicate ratio with tracks (nDuplicateTracks/nAllTracks) = "
          << duplicationRatio_tracks);
  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::INFO,
      "Efficiency with particles (nMatchedParticles/nTrueParticles) = "
          << eff_particle);
  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::INFO,
      "Fake ratio with particles (nFakeParticles/nTrueParticles) = "
          << fakeRatio_particle);
  ACTS_LOG_WITH_LOGGER(
      log, Acts::Logging::INFO,
      "Duplicate ratio with particles (nDuplicateParticles/nTrueParticles) = "
          << duplicationRatio_particle);
}

}  // namespace ActsExamples
