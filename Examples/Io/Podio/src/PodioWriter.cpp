// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioWriter.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"

#include <algorithm>
#include <filesystem>
#include <list>
#include <mutex>

#include <podio/CollectionBase.h>
#include <podio/Frame.h>
#include <tbb/enumerable_thread_specific.h>

namespace ActsExamples {
namespace detail {

using CollectionHandle =
    ConsumeDataHandle<std::unique_ptr<podio::CollectionBase>>;

class PodioWriterImpl {
 public:
  PodioWriterImpl(const PodioWriter::Config& config, PodioWriter& parent)
      : m_cfg(config),
        m_inputPodioFrame(&parent, "InputPodioFrame"),
        m_singleWriter(
            config.separateFilesPerThread
                ? nullptr
                : std::make_unique<ActsPlugins::PodioUtil::ROOTWriter>(
                      config.outputPath)),
        m_useThreadLocalWriters(config.separateFilesPerThread) {}

  PodioWriter::Config m_cfg;

  ConsumeDataHandle<podio::Frame> m_inputPodioFrame;

  std::vector<std::unique_ptr<CollectionHandle>> m_collections;

  std::unique_ptr<ActsPlugins::PodioUtil::ROOTWriter> m_singleWriter;
  tbb::enumerable_thread_specific<
      std::unique_ptr<ActsPlugins::PodioUtil::ROOTWriter>>
      m_threadLocalWriters;

  std::mutex m_writeMutex;
  bool m_useThreadLocalWriters{};

  std::string threadLocalFileName(std::size_t threadId) const {
    namespace fs = std::filesystem;
    const fs::path basePath{m_cfg.outputPath};
    const fs::path directory = basePath.parent_path();
    const std::string stem = basePath.stem().string();
    const std::string extension = basePath.extension().string();

    std::string threadStem = stem.empty() ? basePath.filename().string() : stem;
    threadStem += "_thread" + std::to_string(threadId);

    fs::path threadFile = directory / (threadStem + extension);
    return threadFile.string();
  }

  ActsPlugins::PodioUtil::ROOTWriter& writerForThread(std::size_t threadId) {
    if (!m_useThreadLocalWriters) {
      return *m_singleWriter;
    }
    auto& localWriter = m_threadLocalWriters.local();
    if (!localWriter) {
      localWriter = std::make_unique<ActsPlugins::PodioUtil::ROOTWriter>(
          threadLocalFileName(threadId));
    }
    return *localWriter;
  }
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
    ACTS_DEBUG("- " << collection);
    auto up = std::make_unique<detail::CollectionHandle>(this, collection);
    m_impl->m_collections.push_back(std::move(up));

    m_impl->m_collections.back()->initialize(collection);
  }
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

  std::unique_lock<std::mutex> guard;
  if (!m_impl->m_useThreadLocalWriters) {
    guard = std::unique_lock<std::mutex>(m_impl->m_writeMutex);
  }
  for (const auto& handle : m_impl->m_collections) {
    auto collectionPtr = (*handle)(ctx);
    if (!collectionPtr) {
      ACTS_ERROR("PodioWriter::write: collection is not initialized");
      return ProcessCode::ABORT;
    }
    ACTS_VERBOSE("PodioWriter::write: adding collection " << handle->name()
                                                          << " to frame");
    frame.put(std::move(collectionPtr), handle->name());
  }
  auto& writer = m_impl->writerForThread(ctx.threadId);
  writer.writeFrame(frame, m_impl->m_cfg.category);

  return ProcessCode::SUCCESS;
}

ProcessCode PodioWriter::finalize() {
  if (m_impl->m_useThreadLocalWriters) {
    for (auto& tlsWriter : m_impl->m_threadLocalWriters) {
      if (tlsWriter) {
        tlsWriter->finish();
      }
    }
  } else if (m_impl->m_singleWriter) {
    m_impl->m_singleWriter->finish();
  }

  return ProcessCode::SUCCESS;
}

PodioWriter::~PodioWriter() = default;

const PodioWriter::Config& PodioWriter::config() const {
  return m_impl->m_cfg;
}

}  // namespace ActsExamples
