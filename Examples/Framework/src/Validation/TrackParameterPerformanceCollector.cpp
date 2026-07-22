// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Validation/TrackParameterPerformanceCollector.hpp"

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsExamples/Validation/TruthParametersOnSurface.hpp"

#include <optional>
#include <utility>
#include <vector>

namespace ActsExamples {

namespace {

/// Extract the configured parameters and covariance from a track state. If no
/// parameter type is configured, the best available parameters are used.
std::optional<std::pair<Acts::BoundVector, Acts::BoundMatrix>> stateParameters(
    std::optional<TrackParameterPerformanceCollector::ParameterType>
        parameterType,
    const auto& state) {
  using enum TrackParameterPerformanceCollector::ParameterType;

  if (!parameterType.has_value()) {
    if (!state.hasSmoothed() && !state.hasFiltered() && !state.hasPredicted()) {
      return std::nullopt;
    }
    // best available parameters, i.e. smoothed, filtered, or predicted
    return std::pair(state.parameters(), state.covariance());
  }
  if (parameterType == Predicted && state.hasPredicted()) {
    return std::pair(state.predicted(), state.predictedCovariance());
  }
  if (parameterType == Filtered && state.hasFiltered()) {
    return std::pair(state.filtered(), state.filteredCovariance());
  }
  if (parameterType == Smoothed && state.hasSmoothed()) {
    return std::pair(state.smoothed(), state.smoothedCovariance());
  }
  if (parameterType == Unbiased && state.hasSmoothed() &&
      state.hasProjector() && state.hasCalibrated()) {
    return Acts::calculateUnbiasedParametersCovariance(state);
  }
  return std::nullopt;
}

}  // namespace

ResPlotTool::Config
TrackParameterPerformanceCollector::defaultResPlotToolConfig() {
  using BoostRegularAxis = ResPlotTool::BoostRegularAxis;

  ResPlotTool::Config cfg;
  cfg.paramNames = {"loc0", "loc1", "phi", "theta", "qop", "t"};

  auto& binning = cfg.varBinning;
  binning.erase("Residual_d0");
  binning.erase("Residual_z0");
  binning.insert_or_assign("Residual_loc0",
                           BoostRegularAxis(100, -0.5, 0.5, "r_{loc0} [mm]"));
  binning.insert_or_assign("Residual_loc1",
                           BoostRegularAxis(100, -0.5, 0.5, "r_{loc1} [mm]"));
  // seed-estimated parameters are much coarser than fitted ones
  binning.insert_or_assign("Residual_phi",
                           BoostRegularAxis(100, -0.1, 0.1, "r_{#phi} [rad]"));
  binning.insert_or_assign(
      "Residual_theta", BoostRegularAxis(100, -0.1, 0.1, "r_{#theta} [rad]"));
  binning.insert_or_assign("Residual_qop",
                           BoostRegularAxis(100, -0.5, 0.5, "r_{q/p} [c/GeV]"));
  binning.insert_or_assign(
      "Residual_qopt", BoostRegularAxis(100, -0.5, 0.5, "r_{q/pT} [c/GeV]"));
  binning.insert_or_assign("Residual_qopt_rel",
                           BoostRegularAxis(100, -0.5, 0.5, "r_{rel q/pT}"));

  return cfg;
}

TrackParameterPerformanceCollector::TrackParameterPerformanceCollector(
    Config cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(std::move(cfg)),
      m_logger(std::move(logger)),
      m_resPlotTool(m_cfg.resPlotToolConfig, m_logger->level()) {
  std::vector<Acts::GeometryHierarchyMap<unsigned int>::InputElement> elements;
  elements.reserve(m_cfg.geometrySelection.size());
  for (const Acts::GeometryIdentifier& geoId : m_cfg.geometrySelection) {
    elements.emplace_back(geoId, 0u);
  }
  m_geometrySelection =
      Acts::GeometryHierarchyMap<unsigned int>(std::move(elements));
}

void TrackParameterPerformanceCollector::fill(
    const Acts::GeometryContext& geoContext, const ConstTrackContainer& tracks,
    const SimParticleContainer& particles, const SimHitContainer& simHits,
    const MeasurementParticlesMap& measurementParticlesMap,
    const MeasurementSimHitsMap& measurementSimHitsMap) {
  std::vector<ParticleHitCount> particleHitCounts;

  for (const auto& track : tracks) {
    ++m_stats.nTotalTracks;

    // Require the track to stem from a single truth particle
    identifyContributingParticles(measurementParticlesMap, track,
                                  particleHitCounts);
    if (particleHitCounts.size() != 1) {
      ACTS_VERBOSE("No unique truth particle for track " << track.index());
      continue;
    }

    const auto ip = particles.find(particleHitCounts.front().particleId);
    if (ip == particles.end()) {
      ACTS_DEBUG("Truth particle not found for track " << track.index());
      continue;
    }
    const auto& particle = *ip;

    ++m_stats.nMatchedTracks;

    for (const auto& state : track.trackStatesReversed()) {
      if (!state.typeFlags().isMeasurement() || state.typeFlags().isOutlier()) {
        continue;
      }
      if (!state.hasReferenceSurface()) {
        continue;
      }
      const Acts::Surface& surface = state.referenceSurface();

      if (!m_cfg.geometrySelection.empty() &&
          m_geometrySelection.find(surface.geometryId()) ==
              m_geometrySelection.end()) {
        continue;
      }

      const auto stateParams = stateParameters(m_cfg.parameterType, state);
      if (!stateParams.has_value()) {
        continue;
      }

      const auto* sourceLink =
          state.getUncalibratedSourceLink().getPtr<IndexSourceLink>();
      if (sourceLink == nullptr) {
        continue;
      }

      const std::optional<Acts::BoundTrackParameters> truth =
          truthParametersOnSurface(geoContext, surface, sourceLink->index(),
                                   particle, simHits, measurementSimHitsMap,
                                   logger());
      if (!truth.has_value()) {
        continue;
      }

      ++m_stats.nFilledStates;

      const Acts::BoundTrackParameters parameters(
          surface.getSharedPtr(), stateParams->first, stateParams->second,
          track.particleHypothesis());
      m_resPlotTool.fill(truth.value(), parameters);
    }
  }
}

void TrackParameterPerformanceCollector::logSummary() const {
  ACTS_INFO("=== Track Parameter Performance Summary ===");
  ACTS_INFO("Total tracks: " << m_stats.nTotalTracks);
  ACTS_INFO("Total truth-matched tracks: " << m_stats.nMatchedTracks);
  ACTS_INFO("Total filled track states: " << m_stats.nFilledStates);

  if (m_stats.nTotalTracks > 0) {
    const double fraction =
        static_cast<double>(m_stats.nMatchedTracks) / m_stats.nTotalTracks;
    ACTS_INFO("Truth-matched fraction: " << fraction * 100 << "%");
  }
}

}  // namespace ActsExamples
