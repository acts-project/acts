// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepWriter.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepOutputConverter.hpp"

#include <podio/CollectionBase.h>

namespace ActsExamples {

EDM4hepWriter::EDM4hepWriter(const Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepWriter", level)),
      m_inputPodioFrame(this, "InputPodioFrame"),
      m_writer(config.outputPath) {
  ACTS_VERBOSE("Created output file " << config.outputPath);

  m_inputPodioFrame.initialize(config.inputPodioFrame);

  for (auto& converter : config.converters) {
    converter->initialize(*this);
  }
}

std::string EDM4hepWriter::name() const {
  return "EDM4hepWriter";
}

ProcessCode EDM4hepWriter::write(const AlgorithmContext& ctx) {
  podio::Frame frame = [&]() -> podio::Frame {
    // If we are configured to read a podio::Frame from the event store, do so
    if (m_inputPodioFrame.isInitialized()) {
      return m_inputPodioFrame(ctx);
    } else {
      return podio::Frame();
    }
  }();

  // Invoke converters
  for (auto& converter : m_cfg.converters) {
    converter->convert(ctx, frame);
  }

  std::lock_guard guard(m_writeMutex);
  m_writer.writeFrame(frame, "events");

  return ProcessCode::SUCCESS;
}

ProcessCode EDM4hepWriter::finalize() {
  m_writer.finish();

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
