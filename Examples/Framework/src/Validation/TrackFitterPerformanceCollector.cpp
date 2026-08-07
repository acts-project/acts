// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackFitterPerformanceCollector.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <utility>

namespace ActsExamples {

TrackFitterPerformanceCollector::TrackFitterPerformanceCollector(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(std::move(cfg)),
      m_logger(std::move(logger)),
      m_resPlotTool(m_cfg.resPlotToolConfig, m_logger->level()),
      m_effPlotTool(m_cfg.effPlotToolConfig, m_logger->level()),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig,
                             m_logger->level()) {}

void TrackFitterPerformanceCollector::fill(
    const Acts::GeometryContext& geoContext, const ConstTrackContainer& tracks,
    const SimParticleContainer& particles,
    const TrackParticleMatching& trackParticleMatching) {
  // Truth particles with corresponding reconstructed tracks
  std::vector<SimBarcode> reconParticleIds;
  reconParticleIds.reserve(tracks.size());

  // Loop over all tracks
  for (const auto& track : tracks) {
    ++m_stats.nTotalTracks;

    // Select reco track with fitted parameters
    if (!track.hasReferenceSurface()) {
      ACTS_DEBUG("No fitted track parameters for track " << track.index());
      continue;
    }
    Acts::BoundTrackParameters fittedParameters =
        track.createParametersAtReference();

    // Get the truth-matched particle
    auto imatched = trackParticleMatching.find(track.index());
    if (imatched == trackParticleMatching.end()) {
      ACTS_DEBUG("No truth particle associated with track " << track.index());
      continue;
    }
    const auto& particleMatch = imatched->second;

    if (!particleMatch.particle.has_value()) {
      ACTS_DEBUG("No truth particle associated with track " << track.index());
      continue;
    }

    // Get the barcode of the majority truth particle
    SimBarcode majorityParticleId = particleMatch.particle.value();

    // Find the truth particle via the barcode
    auto ip = particles.find(majorityParticleId);
    if (ip == particles.end()) {
      ACTS_DEBUG("Majority particle not found for track " << track.index());
      continue;
    }

    // Record this majority particle ID
    reconParticleIds.push_back(ip->particleId());

    // Fill residual plots
    m_resPlotTool.fill(geoContext, ip->initialState(), fittedParameters);

    // Fill track summary info
    m_trackSummaryPlotTool.fill(fittedParameters, track.nTrackStates(),
                                track.nMeasurements(), track.nOutliers(),
                                track.nHoles(), track.nSharedHits());
  }

  // Fill the efficiency
  for (const auto& particle : particles) {
    ++m_stats.nTotalParticles;

    bool isReconstructed = false;
    if (Acts::rangeContainsValue(reconParticleIds, particle.particleId())) {
      isReconstructed = true;
      ++m_stats.nTotalMatchedTracks;
      ++m_stats.nTotalMatchedParticles;
    }

    double minDeltaR = -1;
    for (const auto& closeParticle : particles) {
      if (closeParticle.particleId() == particle.particleId()) {
        continue;
      }
      double distance = Acts::VectorHelpers::deltaR(particle.direction(),
                                                    closeParticle.direction());
      if (minDeltaR == -1 || distance < minDeltaR) {
        minDeltaR = distance;
      }
    }

    m_effPlotTool.fill(geoContext, particle.initialState(), minDeltaR,
                       isReconstructed);
  }
}

void TrackFitterPerformanceCollector::logSummary() const {
  ACTS_INFO("=== Track Fitter Performance Summary ===");
  ACTS_INFO("Total tracks: " << m_stats.nTotalTracks);
  ACTS_INFO("Total matched tracks: " << m_stats.nTotalMatchedTracks);
  ACTS_INFO("Total particles: " << m_stats.nTotalParticles);
  ACTS_INFO("Total matched particles: " << m_stats.nTotalMatchedParticles);

  if (m_stats.nTotalTracks > 0) {
    double efficiency =
        static_cast<double>(m_stats.nTotalMatchedTracks) / m_stats.nTotalTracks;
    ACTS_INFO("Track efficiency: " << efficiency * 100 << "%");
  }
}

template <std::size_t Dim>
void TrackFitterPerformanceCollector::addFittedProfiles(
    const std::map<std::string, Acts::Experimental::Histogram<Dim>>& histMap,
    const std::string& meanPrefix, const std::string& widthPrefix,
    std::vector<Acts::Experimental::ValueHistogram<Dim - 1>>& out) const {
  for (const auto& [name, hist] : histMap) {
    // Extract the suffix from the histogram name (e.g., "_d0_vs_eta")
    const std::string& baseName = hist.name();
    const std::string suffix = baseName.substr(baseName.find('_'));

    auto profiles = extractMeanWidthProfiles(
        m_cfg.fitFunction, hist, meanPrefix + suffix, widthPrefix + suffix,
        m_cfg.fitMinEntries, m_cfg.fitSigmaRange, m_cfg.fitIterations,
        logger());
    if (profiles.fitFailureFraction >=
        m_cfg.warningThresholdFitFailureFraction) {
      ACTS_WARNING("Fit failures for " << baseName << ": "
                                       << profiles.fitFailureFraction * 100
                                       << "%");
    }

    out.push_back(std::move(profiles.mean));
    out.push_back(std::move(profiles.width));
  }
}

TrackFitterPerformanceCollector::FittedProfiles
TrackFitterPerformanceCollector::fitProfiles() const {
  FittedProfiles profiles;

  if (!m_cfg.fitFunction) {
    ACTS_WARNING(
        "No fit function configured; skipping mean/width profile "
        "extraction");
    return profiles;
  }

  addFittedProfiles<2>(m_resPlotTool.resVsEta(), "resmean", "reswidth",
                       profiles.profiles1);
  addFittedProfiles<2>(m_resPlotTool.resVsPt(), "resmean", "reswidth",
                       profiles.profiles1);
  addFittedProfiles<3>(m_resPlotTool.resVsEtaPhi(), "resmean", "reswidth",
                       profiles.profiles2);
  addFittedProfiles<3>(m_resPlotTool.resVsEtaPt(), "resmean", "reswidth",
                       profiles.profiles2);

  addFittedProfiles<2>(m_resPlotTool.pullVsEta(), "pullmean", "pullwidth",
                       profiles.profiles1);
  addFittedProfiles<2>(m_resPlotTool.pullVsPt(), "pullmean", "pullwidth",
                       profiles.profiles1);
  addFittedProfiles<3>(m_resPlotTool.pullVsEtaPhi(), "pullmean", "pullwidth",
                       profiles.profiles2);
  addFittedProfiles<3>(m_resPlotTool.pullVsEtaPt(), "pullmean", "pullwidth",
                       profiles.profiles2);

  return profiles;
}

}  // namespace ActsExamples
