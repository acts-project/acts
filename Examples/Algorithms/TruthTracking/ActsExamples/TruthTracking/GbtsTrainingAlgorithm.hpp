// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding2/GbtsTrainingTool.hpp"
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

class GbtsTrainingAlgorithm final : public IAlgorithm {
 public:
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

    Acts::Experimental::GbtsLayerConnectionTool::Config
        gbtsLayerConnectionToolConfig;

    // TO DO: inplace them explicitly into config
    std::string geometryFileDir =
        "/home/ppd/jasper/QT/acts/gbts_layer_geometry.txt";

    std::string outputFileDir =
        "/home/ppd/jasper/QT/acts/connection-tables/"
        "conditional_prob_2kttbar200_new.txt";
  };

  explicit GbtsTrainingAlgorithm(
      const Config& config,
      std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

  ~GbtsTrainingAlgorithm();

 private:
  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};

  ReadDataHandle<ParticleMeasurementsMap> m_inputParticleMeasurementsMap{
      this, "InputParticleMeasurementsMap"};

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};

  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputHits"};

  ReadDataHandle<MeasurementSimHitsMap> m_inputMeasurementSimHitsMap{
      this, "MeasurementSimHitsMap"};

  // make the mutable class thread safe
  mutable std::mutex m_gbtsTrainingToolMutex;
  mutable std::optional<Acts::Experimental::GbtsLayerConnectionTool>
      m_layerConnectionTool;
};

}  // namespace ActsExamples
