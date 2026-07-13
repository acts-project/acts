// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsTrainingTool.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

/// gbst training algorithm for all training or tuning procedures
class GbtsTrainingAlgorithm final : public IAlgorithm {
 public:
  /// Configuration for GBTS training algorithm
  struct Config {
    /// The input truth particles that should be used to create proto tracks.
    std::string inputParticles;
    /// The input particle-measurements map collection.
    std::string inputParticleMeasurementsMap;
    /// The input measurements collection that is used to sort the proto tracks.
    std::string inputMeasurements;
    /// The input sim hits collection that is used to create the proto tracks.
    std::string inputSimHits;
    /// The input measurement-sim hits map collection.
    std::string inputMeasurementSimHitsMap;
    /// The layer connection tool config
    Acts::Experimental::GbtsLayerConnectionTool::Config
        gbtsLayerConnectionToolConfig;
    /// geometry file used for creating layers
    std::string geometryFileDir{};
    /// output directory for layer connection table
    std::string outputFileDir{};
    /// Use legacy athena style table formatting
    bool useOldFormatting = false;
  };

  /// Constructor for GBTS training algorithm
  /// @param config config for GBTS training algorithm
  /// @param logger acts logger
  explicit GbtsTrainingAlgorithm(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// run time execution per event
  /// @param ctx event context
  ProcessCode execute(const AlgorithmContext& ctx) const final;
  /// finalize for tidy up, in this case creates connection table
  ProcessCode finalize() final;
  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// GBTS training algorithm config
  Config m_cfg;
  /// Handle for input truth particles that should be used to create proto
  /// tracks.
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  /// Handle for input particle-measurements map collection
  ReadDataHandle<ParticleMeasurementsMap> m_inputParticleMeasurementsMap{
      this, "InputParticleMeasurementsMap"};
  /// Handle for input measurements collection that is used to sort the proto
  /// tracks
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  /// Handle for input sim hits collection that is used to create the proto
  /// tracks
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputHits"};
  /// Handle for input measurement-sim hits map collection
  ReadDataHandle<MeasurementSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "MeasurementSimHitsMap"};
  /// mutex used for thread safety
  mutable std::mutex m_gbtsTrainingToolMutex;
  /// layer connection tool from core
  mutable std::optional<Acts::Experimental::GbtsLayerConnectionTool>
      m_layerConnectionTool;
};

}  // namespace ActsExamples
