// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IWriter.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <string>

#include <nlohmann/json.hpp>

namespace ActsExamples {

class JsonTrackWriter : public IWriter {
 public:
  struct Config {
    std::string inputTracks;
    std::string inputMeasurementParticlesMap;
    std::string inputClusters;
    std::string inputTrackParticleMatching;

    std::string outputStem;
    std::string outputDir;
  };

  /// Construct the track writer
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  JsonTrackWriter(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~JsonTrackWriter() override;

  std::string name() const override;

  /// Write geometry using the per-event context (optional).
  ProcessCode write(const AlgorithmContext& ctx) override;

  /// Write geometry using the default context.
  ProcessCode finalize() override;

  /// Readonly access to config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
      this, "InputTrackParticleMatching"};
  ReadDataHandle<ClusterContainer> m_inputClusters{this, "InputClusters"};

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
