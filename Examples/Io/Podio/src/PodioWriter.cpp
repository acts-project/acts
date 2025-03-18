// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioWriter.hpp"

#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include <podio/CollectionBase.h>
#include <podio/Frame.h>

namespace ActsExamples {
namespace detail {
class PodioWriterImpl {
 public:
  PodioWriterImpl(const PodioWriter::Config& config, PodioWriter& parent)
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

  PodioWriter::Config m_cfg;

  ConsumeDataHandle<podio::Frame> m_inputPodioFrame;

  std::mutex m_writeMutex;
  Acts::PodioUtil::ROOTWriter m_writer;
};
}  // namespace detail

PodioWriter::PodioWriter(const Config& config, Acts::Logging::Level level)
    : m_logger(Acts::getDefaultLogger("PodioWriter", level)),
      m_impl(std::make_unique<detail::PodioWriterImpl>(config, *this)) {
  ACTS_DEBUG("Created output file " << config.outputPath);
}

std::string PodioWriter::name() const {
  return "PodioWriter";
}

ProcessCode PodioWriter::write(const AlgorithmContext& ctx) {
  podio::Frame frame = m_impl->m_inputPodioFrame(ctx);

  std::lock_guard guard(m_impl->m_writeMutex);
  m_impl->m_writer.writeFrame(frame, m_impl->m_cfg.category);

  return ProcessCode::SUCCESS;
}

ProcessCode PodioWriter::finalize() {
  m_impl->m_writer.finish();

  return ProcessCode::SUCCESS;
}

PodioWriter::~PodioWriter() = default;

const PodioWriter::Config& PodioWriter::config() const {
  return m_impl->m_cfg;
}

}  // namespace ActsExamples
