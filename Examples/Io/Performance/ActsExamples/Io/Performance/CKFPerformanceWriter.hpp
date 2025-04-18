// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/FakeRatePlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"

#include <cstddef>
#include <mutex>
#include <string>

class TFile;
class TTree;
namespace ActsFatras {
class Barcode;
}  // namespace ActsFatras

namespace ActsExamples {
struct AlgorithmContext;

/// Write out the performance of CombinatorialKalmanFilter (CKF), e.g. track
/// efficiency, fake rate etc.
///
/// @TODO: add duplication plots
///
/// A common file can be provided for the writer to attach his TTree, this is
/// done by setting the Config::rootFile pointer to an existing file.
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class CKFPerformanceWriter final : public WriterT<ConstTrackContainer> {
 public:
  struct Config {
    /// Input (found) tracks collection.
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Input track-particle matching.
    std::string inputParticleTrackMatching;
    /// Output filename.
    std::string filePath = "performance_ckf.root";
    /// Output filemode
    std::string fileMode = "RECREATE";

    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
    FakeRatePlotTool::Config fakeRatePlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;

    /// Write additional matching details to a TTree
    bool writeMatchingDetails = false;
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
                     const ConstTrackContainer& tracks) override;

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

  /// For optional output of the matching details
  TTree* m_matchingTree{nullptr};

  /// Variables to fill in the TTree
  uint32_t m_treeEventNr{};
  uint64_t m_treeParticleId{};
  bool m_treeIsMatched{};

  // Adding numbers for efficiency, fake, duplicate calculations
  std::size_t m_nTotalTracks = 0;
  std::size_t m_nTotalMatchedTracks = 0;
  std::size_t m_nTotalFakeTracks = 0;
  std::size_t m_nTotalDuplicateTracks = 0;
  std::size_t m_nTotalParticles = 0;
  std::size_t m_nTotalMatchedParticles = 0;
  std::size_t m_nTotalDuplicateParticles = 0;
  std::size_t m_nTotalFakeParticles = 0;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<ParticleTrackMatching> m_inputParticleTrackMatching{
      this, "InputParticleTrackMatching"};
};

}  // namespace ActsExamples
