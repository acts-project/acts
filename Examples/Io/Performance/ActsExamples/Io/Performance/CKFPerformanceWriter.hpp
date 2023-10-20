// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/FakeRatePlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <cstddef>
#include <functional>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;
namespace ActsFatras {
class Barcode;
}  // namespace ActsFatras

namespace ActsExamples {
struct AlgorithmContext;

/// Write out the performance of CombinatorialKalmanFilter (CKF), e.g.
/// track efficiency, fake rate etc.
/// @TODO: add duplication plots
///
/// A common file can be provided for the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class CKFPerformanceWriter final : public WriterT<TrajectoriesContainer> {
 public:
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;

  struct Config {
    /// Input (found) trajectories collection.
    std::string inputTrajectories;
    /// Input particles collection.
    std::string inputParticles;
    /// Input hit-particles map collection.
    std::string inputMeasurementParticlesMap;
    /// Output filename.
    std::string filePath = "performance_ckf.root";
    /// Output filemode
    std::string fileMode = "RECREATE";
    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
    FakeRatePlotTool::Config fakeRatePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;

    /// Whether to do double matching or not
    bool doubleMatching = false;

    /// Min reco-truth matching probability
    double truthMatchProbMin = 0.5;

    /// function to check if neural network predicted track label is duplicate
    std::function<bool(std::vector<float>&)> duplicatedPredictor = nullptr;
  };

  /// Construct from configuration and log level.
  CKFPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~CKFPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrajectoriesContainer& trajectories) override;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for efficiency
  EffPlotTool m_effPlotTool;
  EffPlotTool::EffPlotCache m_effPlotCache;
  /// Plot tool for fake rate
  FakeRatePlotTool m_fakeRatePlotTool;
  FakeRatePlotTool::FakeRatePlotCache m_fakeRatePlotCache{};
  /// Plot tool for duplication rate
  DuplicationPlotTool m_duplicationPlotTool;
  DuplicationPlotTool::DuplicationPlotCache m_duplicationPlotCache{};
  /// Plot tool for track hit info
  TrackSummaryPlotTool m_trackSummaryPlotTool;
  TrackSummaryPlotTool::TrackSummaryPlotCache m_trackSummaryPlotCache{};

  // Adding numbers for efficiency, fake, duplicate calculations
  size_t m_nTotalTracks = 0;
  size_t m_nTotalMatchedTracks = 0;
  size_t m_nTotalFakeTracks = 0;
  size_t m_nTotalDuplicateTracks = 0;
  size_t m_nTotalParticles = 0;
  size_t m_nTotalMatchedParticles = 0;
  size_t m_nTotalDuplicateParticles = 0;
  size_t m_nTotalFakeParticles = 0;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
};

}  // namespace ActsExamples
