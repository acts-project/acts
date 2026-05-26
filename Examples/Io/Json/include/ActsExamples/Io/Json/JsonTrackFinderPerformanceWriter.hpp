// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/FakePlotTool.hpp"
#include "ActsExamples/Validation/TrackFinderPerformanceCollector.hpp"
#include "ActsExamples/Validation/TrackQualityPlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <map>
#include <mutex>
#include <set>
#include <string>

namespace ActsExamples {

/// Write out the performance of track finding or seeding as a JSON file.
///
/// Produces the same data as RootTrackFinderPerformanceWriter but serialised
/// to JSON using the Acts JSON plugin.  Can therefore be used without ROOT.
/// The same writer covers both CKF track-finding and seeding performance
/// (the latter via addSeedPerformanceWriters) — only the input collections
/// differ.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class JsonTrackFinderPerformanceWriter final
    : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input (found) tracks collection.
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Input particle-track matching.
    std::string inputParticleTrackMatching;
    /// Input particle measurements map.
    std::string inputParticleMeasurementsMap;
    /// Output filename.
    std::string filePath = "performance_ckf.json";

    /// Plot tool configurations (identical to RootTrackFinderPerformanceWriter).
    EffPlotTool::Config effPlotToolConfig;
    FakePlotTool::Config fakePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
    TrackQualityPlotTool::Config trackQualityPlotToolConfig;

    /// Per-subdetector track summary plots, keyed by name.
    std::map<std::string, std::set<int>> subDetectorTrackSummaryVolumes;

    /// Write per-particle matching details as a JSON array.
    bool writeMatchingDetails = false;
  };

  /// Construct from configuration and log level.
  JsonTrackFinderPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~JsonTrackFinderPerformanceWriter() override;

  /// Serialize collected histograms to JSON and write the output file.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters.
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override;

  Config m_cfg;
  std::mutex m_writeMutex;
  TrackFinderPerformanceCollector m_collector;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<ParticleTrackMatching> m_inputParticleTrackMatching{
      this, "InputParticleTrackMatching"};
  ReadDataHandle<ParticleMeasurementsMap> m_inputParticleMeasurementsMap{
      this, "InputParticleMeasurementsMap"};

  /// Per-particle matching records accumulated across all events (optional).
  struct MatchingRecord {
    std::uint32_t eventNr;
    std::uint32_t vertexPrimary;
    std::uint32_t vertexSecondary;
    std::uint32_t particle;
    std::uint32_t generation;
    std::uint32_t subParticle;
    bool isMatched;
  };
  std::vector<MatchingRecord> m_matchingDetails;
};

}  // namespace ActsExamples
