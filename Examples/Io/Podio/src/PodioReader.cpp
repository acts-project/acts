// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioReader.hpp"

#include "Acts/Plugins/Podio/PodioUtil.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include <filesystem>

#include <podio/Frame.h>
#include <tbb/enumerable_thread_specific.h>

namespace ActsExamples {

namespace detail {

class PodioReaderImpl {
 public:
  explicit PodioReaderImpl(PodioReader::Config cfg, PodioReader& parent)
      : m_frameWriteHandle(&parent, "EDM4hepFrameOutput"),
        m_cfg(std::move(cfg)) {
    m_eventsRange = std::make_pair(0, reader().getEntries("events"));
    if (!std::filesystem::exists(m_cfg.inputPath)) {
      throw std::invalid_argument("Input file does not exist");
    }
    if (m_cfg.outputFrame.empty()) {
      throw std::invalid_argument("Output frame name is not set");
    }
    if (m_cfg.category.empty()) {
      throw std::invalid_argument("Category name is not set");
    }

    m_frameWriteHandle.initialize(m_cfg.outputFrame);
  }

  Acts::PodioUtil::ROOTReader& reader() {
    bool exists = false;
    auto& reader = m_reader.local(exists);
    if (!exists) {
      reader.openFile(m_cfg.inputPath);
    }

    return reader;
  }

  WriteDataHandle<podio::Frame> m_frameWriteHandle;
  std::pair<std::size_t, std::size_t> m_eventsRange;

  tbb::enumerable_thread_specific<Acts::PodioUtil::ROOTReader> m_reader;
  PodioReader::Config m_cfg;
};

}  // namespace detail

PodioReader::PodioReader(const Config& config, Acts::Logging::Level level)
    : m_impl(std::make_unique<detail::PodioReaderImpl>(config, *this)),
      m_logger(Acts::getDefaultLogger("PodioReader", level)) {}

PodioReader::~PodioReader() = default;

ProcessCode PodioReader::read(const AlgorithmContext& context) {
  ACTS_DEBUG("Reading EDM4hep inputs");
  podio::Frame frame = m_impl->reader().readEntry(
      m_impl->m_cfg.category, static_cast<unsigned int>(context.eventNumber));
  m_impl->m_frameWriteHandle(context, std::move(frame));

  return ProcessCode::SUCCESS;
}

std::pair<std::size_t, std::size_t> PodioReader::availableEvents() const {
  return m_impl->m_eventsRange;
}

std::string PodioReader::name() const {
  return "PodioReader";
}

const PodioReader::Config& PodioReader::config() const {
  return m_impl->m_cfg;
}

}  // namespace ActsExamples
