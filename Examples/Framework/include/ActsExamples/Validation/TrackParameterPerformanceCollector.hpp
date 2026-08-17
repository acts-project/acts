// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/ParametersOnSurface.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace ActsExamples {

/// Where to take the reconstructed parameters from.
enum class TrackParameterSource {
  /// The track parameters at the track reference surface, i.e. the fitter
  /// output as delivered.
  Track,
  /// The parameters of the individual track states, on the surface they sit
  /// on. Compared against the simulated hits of the state's measurement.
  TrackState,
};

/// Collects performance histograms of the track parameters, without any file
/// I/O.
///
/// Collects residual/pull histograms, efficiency plots, and track summary
/// information for track fitting performance evaluation.
///
/// With `parameterSource = Track` the track parameters at the track reference
/// surface are compared to the truth particle. With `TrackState` every
/// selected measurement state is compared to the truth on its own surface,
/// which is what makes per-sensor estimates, e.g. from a seed, measurable.
///
/// @note The caller must ensure exclusive access (e.g. hold a mutex) when
///       calling fill(). This class applies no locking of its own.
class TrackParameterPerformanceCollector {
 public:
  struct Config {
    ResPlotTool::Config resPlotToolConfig;
    EffPlotTool::Config effPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;

    /// Where to take the reconstructed parameters from.
    TrackParameterSource parameterSource = TrackParameterSource::Track;
    /// Which track-state parameters to use. If not set, the best available
    /// ones (smoothed, filtered, or predicted). `TrackState` source only.
    std::optional<TrackParameterType> parameterType;
    /// If non-empty, only track states in these geometry regions are used.
    /// `TrackState` source only.
    std::vector<Acts::GeometryIdentifier> geometrySelection;

    /// Minimum number of entries in a bin for it to be included in the
    /// mean/width fit.
    int fitMinEntries = 10;
    /// The range in sigma for the iterative Gaussian fit
    double fitSigmaRange = 3.0;
    /// The maximum number of iterations for the iterative Gaussian fit
    int fitIterations = 3;
  };

  TrackParameterPerformanceCollector(
      Config cfg, std::unique_ptr<const Acts::Logger> logger);

  /// Fill histograms for one event.
  ///
  /// @param geoContext the geometry context
  /// @param tracks the input tracks
  /// @param particles the truth particles
  /// @param trackParticleMatching the track to particle matching
  /// @param simHits the simulated hits, required for `TrackState`
  /// @param measurementSimHitsMap the measurement to simulated hits map,
  ///        required for `TrackState`
  ///
  /// @note The caller must ensure exclusive access (e.g. hold a mutex).
  void fill(const Acts::GeometryContext& geoContext,
            const ConstTrackContainer& tracks,
            const SimParticleContainer& particles,
            const TrackParticleMatching& trackParticleMatching,
            const SimHitContainer* simHits = nullptr,
            const MeasurementSimHitsMap* measurementSimHitsMap = nullptr);

  /// Summary count statistics accumulated across all filled events.
  struct Stats {
    std::size_t nTotalTracks = 0;
    std::size_t nTotalMatchedTracks = 0;
    std::size_t nTotalFakeTracks = 0;
    std::size_t nTotalParticles = 0;
    std::size_t nTotalMatchedParticles = 0;
    /// Track states skipped for lack of the requested parameters.
    std::size_t nMissingStateParameters = 0;
    /// Track states skipped for lack of truth hits.
    std::size_t nMissingStateTruth = 0;
  };

  /// Return accumulated event counts.
  const Stats& stats() const { return m_stats; }

  /// Emit summary statistics via the internal logger.
  void logSummary() const;

  /// @name Accessors for the underlying plot tools
  /// @{
  const ResPlotTool& resPlotTool() const { return m_resPlotTool; }
  const EffPlotTool& effPlotTool() const { return m_effPlotTool; }
  const TrackSummaryPlotTool& trackSummaryPlotTool() const {
    return m_trackSummaryPlotTool;
  }
  /// @}

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  /// Fill the residuals of the selected measurement states of one track
  /// against the truth on their own surfaces.
  void fillTrackStates(const Acts::GeometryContext& geoContext,
                       const ConstTrackProxy& track,
                       const SimParticle& particle,
                       const SimHitContainer& simHits,
                       const MeasurementSimHitsMap& measurementSimHitsMap);

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  ResPlotTool m_resPlotTool;
  EffPlotTool m_effPlotTool;
  TrackSummaryPlotTool m_trackSummaryPlotTool;

  Acts::GeometryHierarchyMap<unsigned int> m_geometrySelection;

  Stats m_stats;
};

}  // namespace ActsExamples
