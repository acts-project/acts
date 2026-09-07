// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"
#include "ActsExamples/Validation/TrackParameterPerformanceCollector.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <mutex>
#include <optional>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// Write out the residual and pull of track parameters and efficiency.
///
/// Efficiency here is the fraction of smoothed tracks compared to all tracks.
///
/// With `parameterSource = TrackState` the residuals are taken from the
/// individual measurement states on their own surfaces instead of the track
/// reference surface, optionally restricted to a geometry region. That needs
/// `inputSimHits` and `inputMeasurementSimHitsMap` for the truth.
///
/// With `reference = Measurement` those states are compared against their own
/// calibrated measurement instead. That only constrains the local parameters,
/// but it needs no truth input and therefore also runs on data.
///
/// A common file can be provided for the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootTrackParameterPerformanceWriter final
    : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input track collection.
    std::string inputTracks;
    /// Input particles collection. `Truth` reference only.
    std::string inputParticles;
    /// Input track-particle matching. `Truth` reference only.
    std::string inputTrackParticleMatching;
    /// Input simulated hits collection. `TrackState` source and `Truth`
    /// reference only.
    std::string inputSimHits;
    /// Input measurement to simulated hits map. `TrackState` source and
    /// `Truth` reference only.
    std::string inputMeasurementSimHitsMap;
    /// Output filename.
    std::string filePath = "performance_track_parameters.root";
    /// Plot tool configurations.
    ResPlotTool::Config resPlotToolConfig;
    EffPlotTool::Config effPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;

    /// Where to take the reconstructed parameters from.
    TrackParameterSource parameterSource = TrackParameterSource::Track;
    /// What to compare the reconstructed parameters against. `Measurement` is
    /// `TrackState` source only and requires an explicit @c parameterType.
    TrackParameterReference reference = TrackParameterReference::Truth;
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
    /// Threshold for warning about fit failure fraction in profile extraction.
    double warningThresholdFitFailureFraction = 0.55;
  };

  /// Construct from configuration and log level.
  /// @param config The configuration
  /// @param level The logger level
  RootTrackParameterPerformanceWriter(Config config,
                                      Acts::Logging::Level level);

  ~RootTrackParameterPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const ConstTrackContainer& tracks) override;

  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<MeasurementSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Collector holding all plot tools and per-event counters.
  TrackParameterPerformanceCollector m_collector;
};

}  // namespace ActsExamples
