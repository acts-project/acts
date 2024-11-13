// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IWriter.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <cstddef>
#include <limits>
#include <memory>
#include <string>

namespace ActsExamples {

class JsonTrackWriter : public IWriter {
 public:
  struct Config {
    std::string inputTracks;
    std::string inputMeasurementParticlesMap;
    std::string fileStemp;
    std::string outputDir;
  };

  /// Construct the track writer
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  JsonTrackWriter(const Config& config, Acts::Logging::Level level);

  std::string name() const override;

  /// Write geometry using the per-event context (optional).
  ProcessCode write(const AlgorithmContext& ctx) override;

  /// Write geometry using the default context.
  ProcessCode finalize() override { return ProcessCode::SUCCESS; }

  /// Readonly access to config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<HitParticlesMap> m_inputMeasurementParticlesMap{
      this, "InputMeasurementParticlesMap"};

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
