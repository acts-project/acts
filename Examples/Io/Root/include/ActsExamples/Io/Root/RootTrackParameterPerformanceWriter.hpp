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
#include "ActsExamples/Validation/ResPlotTool.hpp"
#include "ActsExamples/Validation/TrackParameterPerformanceCollector.hpp"

#include <mutex>
#include <optional>
#include <string>
#include <vector>

class TFile;

namespace ActsExamples {

/// Write out the residual and pull of track parameters against the truth hit
/// information on the reference surface.
///
/// This is intended for track parameters estimated from seeds, e.g. produced
/// by `SeedsToTracks`, which are expressed on the surface of the innermost
/// space point. The input tracks can be selected upstream, e.g. by layer with
/// `TrackSelectorAlgorithm` and geometry-based measurement requirements.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootTrackParameterPerformanceWriter final
    : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input track collection holding the estimated track parameters.
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input simulated hits collection.
    std::string inputSimHits;
    /// Input measurement to particles map.
    std::string inputMeasurementParticlesMap;
    /// Input measurement to simulated hits map.
    std::string inputMeasurementSimHitsMap;
    /// Output filename.
    std::string filePath = "performance_track_parameters.root";
    /// Plot tool configuration.
    ResPlotTool::Config resPlotToolConfig =
        TrackParameterPerformanceCollector::defaultResPlotToolConfig();
    /// The track-state parameters to compare to truth. If not set, the best
    /// available parameters are used (smoothed, filtered, or predicted).
    std::optional<TrackParameterType> parameterType;
    /// Optional geometry selection of track states. If non-empty, only track
    /// states within the given geometry hierarchy regions are used.
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
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  ReadDataHandle<MeasurementSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "InputMeasurementSimHitsMap"};

  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Collector holding the plot tool and per-event counters.
  TrackParameterPerformanceCollector m_collector;
};

}  // namespace ActsExamples
