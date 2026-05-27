// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <mutex>
#include <string>

namespace ActsExamples {

/// Write out residual, pull and efficiency plots for fitted tracks as JSON.
///
/// Mirrors RootTrackFitterPerformanceWriter exactly except:
///  - output is a JSON file (no ROOT dependency)
///  - the Gaussian profile-extraction step (extractMeanWidthProfiles) is
///    omitted; the raw 2D/3D histograms contain strictly more information and
///    profile extraction can be performed downstream in Python.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class JsonTrackFitterPerformanceWriter final
    : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input (fitted) track collection.
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Output filename.
    std::string filePath = "performance_track_fitter.json";

    /// Plot tool configurations (identical to
    /// RootTrackFitterPerformanceWriter).
    ResPlotTool::Config resPlotToolConfig;
    EffPlotTool::Config effPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;

    /// Kept for API compatibility with RootTrackFitterPerformanceWriter;
    /// not used in the JSON writer (no in-C++ Gaussian fitting).
    int fitMinEntries = 10;
    double fitSigmaRange = 3.0;
    int fitIterations = 3;
    double warningThresholdFitFailureFraction = 0.55;
  };

  /// Construct from configuration and log level.
  JsonTrackFitterPerformanceWriter(Config config, Acts::Logging::Level level);

  ~JsonTrackFitterPerformanceWriter() override;

  /// Serialize collected histograms to JSON and write the output file.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters.
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override;

  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};

  std::mutex m_writeMutex;
  ResPlotTool m_resPlotTool;
  EffPlotTool m_effPlotTool;
  TrackSummaryPlotTool m_trackSummaryPlotTool;
};

}  // namespace ActsExamples
