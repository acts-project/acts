// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/DuplicationPlotTool.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"

#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

class SeedingPerformanceWriter final : public WriterT<SimSeedContainer> {
 public:
  using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;

  struct Config {
    /// Input reconstructed seed collection.
    std::string inputSeeds;
    /// Input hit to particles map.
    std::string inputMeasurementParticlesMap;
    /// Input truth particles collection.
    std::string inputParticles;
    /// Output filename.
    std::string filePath = "performance_track_seeding.root";
    /// Output file mode
    std::string fileMode = "RECREATE";
    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
    DuplicationPlotTool::Config duplicationPlotToolConfig;
  };

  /// Construct from configuration and log level.
  /// @param config The configuration
  /// @param level
  SeedingPerformanceWriter(Config config, Acts::Logging::Level level);

  ~SeedingPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SimSeedContainer& seeds) override;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for efficiency
  EffPlotTool m_effPlotTool;
  EffPlotTool::EffPlotCache m_effPlotCache;
  /// Plot tool for duplication rate
  DuplicationPlotTool m_duplicationPlotTool;
  DuplicationPlotTool::DuplicationPlotCache m_duplicationPlotCache{};

  size_t m_nTotalSeeds = 0;
  size_t m_nTotalMatchedSeeds = 0;
  size_t m_nTotalParticles = 0;
  size_t m_nTotalMatchedParticles = 0;
  size_t m_nTotalDuplicatedParticles = 0;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMaps"};
};

}  // namespace ActsExamples
