// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Validation/DuplicationPlotTool.hpp"
#include "ACTFW/Validation/EffPlotTool.hpp"
#include "ACTFW/Validation/FakeRatePlotTool.hpp"
#include "ACTFW/Validation/TrackSummaryPlotTool.hpp"
#include "Acts/Utilities/Units.hpp"

class TFile;
class TTree;

using namespace Acts::UnitLiterals;

namespace FW {

/// Write out the performance of CombinatorialKalmanFilter (CKF), e.g.
/// track efficiency, fake rate etc.
/// @TODO: add duplication plots
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class CKFPerformanceWriter final : public WriterT<TrajectoryContainer> {
 public:
  struct Config {
    /// Input truth particles collection.
    std::string inputParticles;
    /// Input (found) trajectories collection.
    std::string inputTrajectories;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "performance_ckf.root";
    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
    FakeRatePlotTool::Config fakeRatePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
    /// Min reco-truth matching probability
    double truthMatchProbMin = 0.5;
    /// Min number of measurements
    size_t nMeasurementsMin = 9;
    /// Min transverse momentum
    double ptMin = 1_GeV;
  };

  /// Construct from configuration and log level.
  CKFPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~CKFPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode endRun() final override;

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const TrajectoryContainer& trajectories) final override;

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
  TrackSummaryPlotTool::TrackSummaryPlotCache m_trackSummaryPlotCache;
};

}  // namespace FW
