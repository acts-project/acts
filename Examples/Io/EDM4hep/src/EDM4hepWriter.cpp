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
namespace detail {
class EDM4hepWriterImpl {
 public:
  EDM4hepWriterImpl(const EDM4hepWriter::Config& config, EDM4hepWriter& parent)
      : m_cfg(config),
        m_inputPodioFrame(&parent, "InputPodioFrame"),
        m_writer(config.outputPath) {
    if (m_cfg.category.empty()) {
      throw std::invalid_argument("Category name is not set");
    }
    if (m_cfg.inputFrame.empty()) {
      throw std::invalid_argument("Input frame name is not set");
    }
    m_inputPodioFrame.initialize(m_cfg.inputFrame);
  }

  EDM4hepWriter::Config m_cfg;

  ConsumeDataHandle<podio::Frame> m_inputPodioFrame;

  std::mutex m_writeMutex;
  Acts::PodioUtil::ROOTWriter m_writer;
};
}  // namespace detail

EDM4hepWriter::EDM4hepWriter(const Config& config, Acts::Logging::Level level)
    : m_logger(Acts::getDefaultLogger("EDM4hepWriter", level)),
      m_impl(std::make_unique<detail::EDM4hepWriterImpl>(config, *this)) {
  ACTS_DEBUG("Created output file " << config.outputPath);
}

std::string EDM4hepWriter::name() const {
  return "EDM4hepWriter";
}

ProcessCode EDM4hepWriter::write(const AlgorithmContext& ctx) {
  podio::Frame frame = m_impl->m_inputPodioFrame(ctx);

  std::lock_guard guard(m_impl->m_writeMutex);
  m_impl->m_writer.writeFrame(frame, m_impl->m_cfg.category);

  return ProcessCode::SUCCESS;
}

ProcessCode EDM4hepWriter::finalize() {
  m_impl->m_writer.finish();

  return ProcessCode::SUCCESS;
}

EDM4hepWriter::~EDM4hepWriter() = default;

const EDM4hepWriter::Config& EDM4hepWriter::config() const {
  return m_impl->m_cfg;
}

}  // namespace ActsExamples
