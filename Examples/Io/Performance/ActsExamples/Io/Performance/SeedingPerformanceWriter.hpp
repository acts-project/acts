// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Seeding/SimSpacePoint.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/FakeRatePlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <mutex>

class TFile;
class TTree;

// using namespace Acts::UnitLiterals;

namespace ActsExamples {

class SeedingPerformanceWriter final
    : public WriterT<std::vector<std::vector<Acts::Seed<SimSpacePoint>>>> {
 public:
  struct Config {
    /// Input hit to particles map
    std::string inputHitParticlesMap;
    /// Input truth particles collection.
    std::string inputParticles;
    /// Input seeds to be analyzed.
    std::string inputSeeds;
    /// Input seeds as proto tracks
    std::string inputProtoTracks;
    /// input Clusters from the event#-hits.csv file for calculating efficiency.
    std::string inputClusters;
    /// Output filename.
    std::string outputFilename = "performance_track_seeding.root";

    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
    FakeRatePlotTool::Config fakeRatePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
  };

  /// @brief Finds all the particles that are in common to all space points in
  /// the seed.
  /// @param seed The seed to be analyzed
  std::set<ActsFatras::Barcode> identifySharedParticles(
      const Acts::Seed<SimSpacePoint>* seed) const;

  /// Construct from configuration and log level.
  SeedingPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~SeedingPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode endRun() final override;

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::vector<std::vector<Acts::Seed<SimSpacePoint>>>&
                         seedVector) final override;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for efficiency
  EffPlotTool m_effPlotTool;
  EffPlotTool::EffPlotCache m_effPlotCache;
  // /// Plot tool for fake rate
  // FakeRatePlotTool m_fakeRatePlotTool;
  // FakeRatePlotTool::FakeRatePlotCache m_fakeRatePlotCache{};
  // // /// Plot tool for duplication rate
  // DuplicationPlotTool m_duplicationPlotTool;
  // DuplicationPlotTool::DuplicationPlotCache m_duplicationPlotCache{};
  // /// Plot tool for track hit info
  // TrackSummaryPlotTool m_trackSummaryPlotTool;
  // TrackSummaryPlotTool::TrackSummaryPlotCache m_trackSummaryPlotCache;
};

}  // namespace ActsExamples
