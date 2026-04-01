// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/FakePlotTool.hpp"
#include "ActsExamples/Validation/TrackQualityPlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <cstddef>
#include <map>
#include <set>
#include <string>

namespace ActsExamples {

struct AlgorithmContext;

/// Collects track-finder performance histograms without any file I/O.
///
/// Fill histograms for each event via fill(). Caller is responsible for
/// thread safety — this class applies no locking of its own.
class TrackFinderPerformanceCollector {
 public:
  struct Config {
    EffPlotTool::Config effPlotToolConfig;
    FakePlotTool::Config fakePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
    TrackQualityPlotTool::Config trackQualityPlotToolConfig;

    /// Optional per-subdetector track summary plots, keyed by name.
    /// The value is the set of geometry volume IDs to include.
    std::map<std::string, std::set<int>> subDetectorTrackSummaryVolumes;
  };

  TrackFinderPerformanceCollector(Config cfg,
                                  std::unique_ptr<const Acts::Logger> logger);

  /// Fill histograms for one event.
  ///
  /// @note The caller must ensure exclusive access (e.g. hold a mutex).
  ProcessCode fill(const AlgorithmContext& ctx,
                   const ConstTrackContainer& tracks,
                   const SimParticleContainer& particles,
                   const TrackParticleMatching& trackParticleMatching,
                   const ParticleTrackMatching& particleTrackMatching,
                   const InverseMultimap<SimBarcode>& particleMeasurementsMap);

  /// Summary count statistics accumulated across all filled events.
  struct Stats {
    std::size_t nTotalTracks = 0;
    std::size_t nTotalMatchedTracks = 0;
    std::size_t nTotalFakeTracks = 0;
    std::size_t nTotalDuplicateTracks = 0;
    std::size_t nTotalParticles = 0;
    std::size_t nTotalMatchedParticles = 0;
    std::size_t nTotalDuplicateParticles = 0;
    std::size_t nTotalFakeParticles = 0;
  };

  /// Return accumulated event counts.
  const Stats& stats() const { return m_stats; }

  /// Emit efficiency/fake/duplicate summary statistics via the internal logger.
  void logSummary() const;

  /// @name Accessors for the underlying plot tools
  /// @{
  const EffPlotTool& effPlotTool() const { return m_effPlotTool; }
  const FakePlotTool& fakePlotTool() const { return m_fakePlotTool; }
  const DuplicationPlotTool& duplicationPlotTool() const {
    return m_duplicationPlotTool;
  }
  const TrackSummaryPlotTool& trackSummaryPlotTool() const {
    return m_trackSummaryPlotTool;
  }
  const std::map<std::string, TrackSummaryPlotTool>& subDetectorSummaryTools()
      const {
    return m_subDetectorSummaryTools;
  }
  const TrackQualityPlotTool& trackQualityPlotTool() const {
    return m_trackQualityPlotTool;
  }
  /// @}

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  EffPlotTool m_effPlotTool;
  FakePlotTool m_fakePlotTool;
  DuplicationPlotTool m_duplicationPlotTool;
  TrackSummaryPlotTool m_trackSummaryPlotTool;
  std::map<std::string, TrackSummaryPlotTool> m_subDetectorSummaryTools;
  TrackQualityPlotTool m_trackQualityPlotTool;

  Stats m_stats;
};

}  // namespace ActsExamples
