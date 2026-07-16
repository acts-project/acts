// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <cstddef>
#include <map>
#include <string>

namespace ActsExamples {

/// Collects track-fitter performance histograms without any file I/O.
///
/// Collects residual/pull histograms, efficiency plots, and track summary
/// information for track fitting performance evaluation.
///
/// @note The caller must ensure exclusive access (e.g. hold a mutex) when
///       calling fill(). This class applies no locking of its own.
class TrackFitterPerformanceCollector {
 public:
  struct Config {
    ResPlotTool::Config resPlotToolConfig;
    EffPlotTool::Config effPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;

    /// Minimum number of entries in a bin for it to be included in the
    /// mean/width fit.
    int fitMinEntries = 10;
    /// The range in sigma for the iterative Gaussian fit
    double fitSigmaRange = 3.0;
    /// The maximum number of iterations for the iterative Gaussian fit
    int fitIterations = 3;
  };

  TrackFitterPerformanceCollector(Config cfg,
                                  std::unique_ptr<const Acts::Logger> logger);

  /// Fill histograms for one event.
  ///
  /// @note The caller must ensure exclusive access (e.g. hold a mutex).
  void fill(const Acts::GeometryContext& geoContext,
            const ConstTrackContainer& tracks,
            const SimParticleContainer& particles,
            const TrackParticleMatching& trackParticleMatching);

  /// Summary count statistics accumulated across all filled events.
  struct Stats {
    std::size_t nTotalTracks = 0;
    std::size_t nTotalMatchedTracks = 0;
    std::size_t nTotalFakeTracks = 0;
    std::size_t nTotalParticles = 0;
    std::size_t nTotalMatchedParticles = 0;
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

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  ResPlotTool m_resPlotTool;
  EffPlotTool m_effPlotTool;
  TrackSummaryPlotTool m_trackSummaryPlotTool;

  Stats m_stats;
};

}  // namespace ActsExamples
