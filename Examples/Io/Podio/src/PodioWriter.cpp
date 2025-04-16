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

#include <algorithm>
#include <list>

#include <podio/CollectionBase.h>
#include <podio/Frame.h>

namespace ActsExamples {
namespace detail {

using CollectionHandle =
    ConsumeDataHandle<std::unique_ptr<podio::CollectionBase>>;

class PodioWriterImpl {
 public:
  PodioWriterImpl(const PodioWriter::Config& config, PodioWriter& parent)
      : m_cfg(config),
        m_inputPodioFrame(&parent, "InputPodioFrame"),
        m_writer(config.outputPath) {}

  PodioWriter::Config m_cfg;

  ConsumeDataHandle<podio::Frame> m_inputPodioFrame;

  std::vector<std::unique_ptr<CollectionHandle>> m_collections;

  std::mutex m_writeMutex;
  Acts::PodioUtil::ROOTWriter m_writer;
};
}  // namespace detail

PodioWriter::PodioWriter(const Config& config, Acts::Logging::Level level)
    : m_logger(Acts::getDefaultLogger("PodioWriter", level)),
      m_impl(std::make_unique<detail::PodioWriterImpl>(config, *this)) {
  ACTS_DEBUG("Creating output file " << config.outputPath);

  if (m_impl->m_cfg.category.empty()) {
    throw std::invalid_argument("Category name is not set");
  }
  if (!m_impl->m_cfg.inputFrame.has_value()) {
    ACTS_DEBUG("No input frame name set, will create a new one");
  } else {
    m_impl->m_inputPodioFrame.initialize(m_impl->m_cfg.inputFrame.value());
  }

  std::set<std::string> collections;
  std::ranges::copy(m_impl->m_cfg.collections,
                    std::inserter(collections, collections.end()));
  if (collections.size() != m_impl->m_cfg.collections.size()) {
    throw std::invalid_argument(
        "Duplicate collection names in config, please use check your "
        "configuration");
  }
  if (std::ranges::any_of(m_impl->m_cfg.collections,
                          [](const auto& c) { return c.empty(); })) {
    throw std::invalid_argument("Collection name is empty");
  }

  ACTS_DEBUG("Adding " << m_impl->m_cfg.collections.size() << " collections");
  for (const auto& collection : m_impl->m_cfg.collections) {
    auto up = std::make_unique<detail::CollectionHandle>(this, collection);
    m_impl->m_collections.push_back(std::move(up));

    m_impl->m_collections.back()->initialize(collection);
  }
  ACTS_DEBUG("PodioWriter::PodioWriter: " << m_impl->m_collections.size());
}

std::string PodioWriter::name() const {
  return "PodioWriter";
}

ProcessCode PodioWriter::write(const AlgorithmContext& ctx) {
  ACTS_DEBUG("PodioWriter::write");
  podio::Frame frame = [this, &ctx]() -> podio::Frame {
    if (m_impl->m_inputPodioFrame.isInitialized()) {
      ACTS_VERBOSE(
          "PodioWriter::write: taking inputPodioFrame from WhiteBoard");
      return m_impl->m_inputPodioFrame(ctx);
    } else {
      ACTS_VERBOSE("PodioWriter::write: creating new inputPodioFrame");
      return {};
    }
  }();

  std::lock_guard guard(m_impl->m_writeMutex);
  for (const auto& handle : m_impl->m_collections) {
    auto collectionPtr = (*handle)(ctx);
    if (!collectionPtr) {
      ACTS_ERROR("PodioWriter::write: collection is not initialized");
      return ProcessCode::ABORT;
    }
    frame.put(std::move(collectionPtr), handle->name());
  }
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
