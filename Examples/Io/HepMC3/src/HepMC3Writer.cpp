// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"

#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Metadata.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <filesystem>
#include <stdexcept>

#include <HepMC3/Version.h>
#include <HepMC3/WriterAscii.h>
#include <boost/algorithm/string/join.hpp>

namespace ActsExamples {

HepMC3Writer::HepMC3Writer(const Config& config, Acts::Logging::Level level)
    : WriterT(config.inputEvent, "HepMC3Writer", level),
      m_cfg(config),
      m_queueSemaphore{static_cast<long>(m_cfg.maxEventsPending + 1)} {
  if (m_cfg.outputPath.empty()) {
    throw std::invalid_argument("Missing output file path");
  }

  // Resolve compression from config and/or path
  auto [compression, basePath] = resolveCompression(m_cfg.outputPath);
  m_compression = compression;

  // Validate the resolved compression is supported
  if (std::ranges::none_of(HepMC3Util::availableCompressionModes(),
                           [this](HepMC3Util::Compression c) {
                             return c == this->m_compression;
                           })) {
    std::stringstream ss;
    ss << "Unsupported compression mode: " << m_compression;
    throw std::invalid_argument(ss.str());
  }

  // Compute the actual output path with compression extension
  auto format = HepMC3Util::formatFromFilename(basePath);
  if (format == HepMC3Util::Format::root) {
    m_actualOutputPath = basePath;
  } else {
    m_actualOutputPath =
        basePath.string() +
        std::string{HepMC3Util::compressionExtension(m_compression)};
  }

  // Validate output path
  auto absolute = std::filesystem::absolute(m_actualOutputPath);
  if (std::filesystem::exists(absolute) &&
      std::filesystem::is_directory(absolute)) {
    throw std::invalid_argument("Output path is a directory: " +
                                absolute.string());
  }

  if (!std::filesystem::exists(absolute.parent_path())) {
    throw std::invalid_argument("Directory to write into does not exist: " +
                                absolute.parent_path().string());
  }

  // Create the writer
  m_writer =
      HepMC3Util::createWriter(m_actualOutputPath, format, m_compression);
}

HepMC3Writer::~HepMC3Writer() = default;

std::pair<HepMC3Util::Compression, std::filesystem::path>
HepMC3Writer::resolveCompression(const std::filesystem::path& path) const {
  using Compression = HepMC3Util::Compression;

  // Deduce compression from the path
  Compression pathCompression = HepMC3Util::compressionFromFilename(path);

  // Get the base path without compression extension
  std::filesystem::path basePath = path;
  std::string_view ext = HepMC3Util::compressionExtension(pathCompression);
  if (!ext.empty() && path.string().ends_with(ext)) {
    // Remove the compression extension
    basePath = path.string().substr(0, path.string().size() - ext.size());
  }

  // Resolve the compression mode
  Compression resolvedCompression = Compression::none;

  if (m_cfg.compression.has_value()) {
    // Compression explicitly set in config
    resolvedCompression = m_cfg.compression.value();

    // Check consistency if path also specifies compression
    if (pathCompression != Compression::none &&
        pathCompression != resolvedCompression) {
      std::stringstream ss;
      ss << "Compression mismatch: config specifies " << resolvedCompression
         << ", but path '" << path << "' implies " << pathCompression;
      throw std::invalid_argument(ss.str());
    }
  } else {
    // Deduce compression from path
    resolvedCompression = pathCompression;
  }

  return {resolvedCompression, basePath};
}

ProcessCode HepMC3Writer::beginEvent(std::size_t threadId) {
  ACTS_VERBOSE("Begin event, next_event=" << m_nextEvent);
  if (!m_cfg.writeEventsInOrder) {
    // Nothing to do if we don't write in order
    return ProcessCode::SUCCESS;
  }

  ACTS_DEBUG("thread=" << threadId << ", next_event=" << m_nextEvent
                       << " waiting for semaphore");
  m_waiting++;
  std::chrono::seconds timeout{10};
  if (!m_queueSemaphore.try_acquire_for(timeout)) {
    ACTS_ERROR("thread=" << threadId << ", next_event=" << m_nextEvent
                         << " failed to acquire semaphore after "
                         << timeout.count() << "s");
    return ProcessCode::ABORT;
  }
  m_waiting--;
  ACTS_DEBUG("thread=" << threadId << ", next_event=" << m_nextEvent
                       << " have semaphore");

  return ProcessCode::SUCCESS;
}

ProcessCode HepMC3Writer::writeT(
    const AlgorithmContext& ctx,
    const std::shared_ptr<HepMC3::GenEvent>& event) {
  ACTS_VERBOSE("Write: event_nr=" << ctx.eventNumber
                                  << ", thread=" << ctx.threadId);

  if (!m_cfg.writeEventsInOrder) {
    std::scoped_lock lock{m_mutex};
    // Unconditionally write events in whatever order they come in
    m_writer->write_event(*event);
    if (m_writer->failed()) {
      ACTS_ERROR("Failed to write event number: " << ctx.eventNumber);
      return ProcessCode::ABORT;
    }
    m_eventsWritten++;
    return ProcessCode::SUCCESS;
  }

  std::scoped_lock lock{m_mutex};

  auto printQueue = [&]() {
    ACTS_VERBOSE("queue=[" << [&]() {
      std::vector<std::string> numbers;
      numbers.reserve(m_eventQueue.size());
      std::ranges::transform(
          m_eventQueue, std::back_inserter(numbers),
          [](const auto& pair) { return std::to_string(pair.first); });

      return boost::algorithm::join(numbers, ", ");
    }() << "]");
  };

  if (ctx.eventNumber == m_nextEvent) {
    ACTS_DEBUG("event_nr=" << ctx.eventNumber
                           << " is the next event -> writing");

    // write
    m_writer->write_event(*event);
    if (m_writer->failed()) {
      ACTS_ERROR("Failed to write event number: " << ctx.eventNumber);
      return ProcessCode::ABORT;
    }
    m_nextEvent++;

    std::size_t nWritten = 1;

    printQueue();

    while (!m_eventQueue.empty() &&
           m_eventQueue.front().first == static_cast<long long>(m_nextEvent)) {
      auto [nextEventNumber, nextEvent] = std::move(m_eventQueue.front());
      ACTS_VERBOSE("Writing event number: " << nextEventNumber);
      m_eventQueue.erase(m_eventQueue.begin());

      m_writer->write_event(*nextEvent);
      if (m_writer->failed()) {
        ACTS_ERROR("Failed to write event number: " << nextEventNumber);
        return ProcessCode::ABORT;
      }
      m_nextEvent++;
      nWritten++;
    }

    m_eventsWritten += nWritten;

    ACTS_VERBOSE("Wrote " << nWritten << " events, next_event=" << m_nextEvent
                          << ", thread=" << ctx.threadId);
    m_queueSemaphore.release(nWritten);
    ACTS_VERBOSE("thread=" << ctx.threadId << ", released n=" << nWritten
                           << ", waiting=" << m_waiting);

  } else {
    ACTS_DEBUG("event_nr=" << ctx.eventNumber
                           << " is not the next event -> queueing");

    ACTS_VERBOSE(
        "Finding insert location for event number: " << ctx.eventNumber);
    auto it = std::ranges::upper_bound(m_eventQueue, ctx.eventNumber, {},
                                       [](const auto& v) { return v.first; });
    if (it == m_eventQueue.end()) {
      ACTS_VERBOSE("Insert location for " << ctx.eventNumber
                                          << " is at the end of the queue");
    } else {
      ACTS_VERBOSE("Insert location for "
                   << ctx.eventNumber
                   << " is before event number: " << it->first);
    }

    m_eventQueue.emplace(it, ctx.eventNumber, event);
    printQueue();
    m_maxEventQueueSize = std::max(m_maxEventQueueSize, m_eventQueue.size());
  }

  return ProcessCode::SUCCESS;
}

ProcessCode HepMC3Writer::finalize() {
  ACTS_VERBOSE("Finalizing HepMC3Writer");
  if (m_writer) {
    m_writer->close();
  }
  ACTS_DEBUG("max_queue_size=" << m_maxEventQueueSize
                               << " limit=" << m_cfg.maxEventsPending);

  // Write sidecar metadata file with event count
  ACTS_DEBUG("Writing sidecar metadata for " << m_actualOutputPath << " with "
                                             << m_eventsWritten << " events");
  if (!HepMC3Metadata::writeSidecar(
          m_actualOutputPath,
          HepMC3Metadata::HepMC3Metadata{.numEvents = m_eventsWritten},
          logger())) {
    ACTS_WARNING("Failed to write sidecar metadata file for "
                 << m_actualOutputPath);
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
