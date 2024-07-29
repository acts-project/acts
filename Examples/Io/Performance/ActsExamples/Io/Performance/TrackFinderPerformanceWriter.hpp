// This file is part of the Acts project.
//
// Copyright (C) 2019-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/FakeRatePlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <memory>
#include <string>

namespace ActsExamples {
struct AlgorithmContext;

/// Write track finder performance measures.
///
/// Only considers the track finding itself, i.e. grouping of hits into tracks,
/// and computes relevant per-track and per-particles statistics.
class TrackFinderPerformanceWriter final : public WriterT<TrackContainer> {
 public:
  struct Config {
    /// Input reconstructed tracks collection.
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Input particle-track matching.
    std::string inputParticleTrackMatching;
    /// Output filename.
    std::string filePath = "performance_track_finder.root";
    /// Output file mode
    std::string fileMode = "RECREATE";
    /// Output tree name for the tracks
    std::string treeNameTracks = "track_finder_tracks";
    /// Output tree name for the particles
    std::string treeNameParticles = "track_finder_particles";

    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
    FakeRatePlotTool::Config fakeRatePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;

    /// Write additional matching details to a TTree
    bool writeMatchingDetails = false;
  };

  /// Constructor
  /// @param config the configuration
  /// @param level The log level
  TrackFinderPerformanceWriter(const Config& config,
                               Acts::Logging::Level level);

  ~TrackFinderPerformanceWriter() override;

  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const;

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrackContainer& tracks) override;

  struct Impl;

  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples
