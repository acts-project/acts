// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IWriter.hpp"

#include <memory>
#include <string>

#include <DD4hep/DetElement.h>
#include <edm4hep/MCParticleCollection.h>
#include <tbb/enumerable_thread_specific.h>

namespace podio {
class Frame;
}

namespace ActsExamples {

class EDM4hepOutputConverter;

class EDM4hepWriter final : public IWriter {
 public:
  struct Config {
    std::string outputPath;

    /// Retrieve a @c podio::Frame from the event store using this name.
    /// If empty, a new frame will be created by this writer.
    std::string inputPodioFrame = "";

    std::vector<std::shared_ptr<EDM4hepOutputConverter>> converters;
  };

  EDM4hepWriter(const Config& config, Acts::Logging::Level level);

  std::string name() const final;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

  ProcessCode finalize() final;

  ProcessCode write(const AlgorithmContext& ctx) final;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;

  ConsumeDataHandle<podio::Frame> m_inputPodioFrame;

  Acts::PodioUtil::ROOTWriter m_writer;
  std::mutex m_writeMutex;
};

}  // namespace ActsExamples
