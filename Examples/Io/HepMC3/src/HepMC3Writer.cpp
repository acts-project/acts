// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"

#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <filesystem>
#include <stdexcept>

#include <HepMC3/Version.h>
#include <HepMC3/WriterAscii.h>

#if HEPMC3_VERSION_CODE == 3002007
#include "./CompressedIO.h"
#endif

#ifdef HEPMC3_USE_COMPRESSION
#include <HepMC3/WriterGZ.h>
#endif

#ifdef ACTS_HEPMC3_ROOT_SUPPORT
#include <HepMC3/WriterRootTree.h>
#endif

#include <boost/algorithm/string/join.hpp>

namespace ActsExamples {

HepMC3Writer::HepMC3Writer(const Config& config, Acts::Logging::Level level)
    : WriterT(config.inputEvent, "HepMC3Writer", level),
      m_cfg(config),
      m_queueSemaphore{static_cast<long>(m_cfg.maxEventsPending + 1)} {
  if (m_cfg.outputPath.empty()) {
    throw std::invalid_argument("Missing output file path");
  }

  if (std::ranges::none_of(HepMC3Util::availableCompressionModes(),
                           [this](HepMC3Util::Compression c) {
                             return c == this->m_cfg.compression;
                           })) {
    std::stringstream ss;
    ss << "Unsupported compression mode: " << m_cfg.compression;
    throw std::invalid_argument(ss.str());
  }

  if (!m_cfg.perEvent) {
    auto absolute = std::filesystem::absolute(m_cfg.outputPath);
    if (std::filesystem::exists(absolute) &&
        std::filesystem::is_directory(absolute)) {
      throw std::invalid_argument("Output path is a directory: " +
                                  absolute.string());
    }

    if (!std::filesystem::exists(absolute.parent_path())) {
      throw std::invalid_argument("Directory to write into does not exist: " +
                                  absolute.parent_path().string());
    }
    // Create a single file writer
    m_writer = createWriter(m_cfg.outputPath);
  }

  if (m_cfg.perEvent &&
      HepMC3Util::formatFromFilename(m_cfg.outputPath.string()) ==
          HepMC3Util::Format::root) {
    ACTS_WARNING(
        "Per-event output is enabled and the output format is ROOT. This is "
        "likely not what you want");
  }
}

std::unique_ptr<HepMC3::Writer> HepMC3Writer::createWriter(
    const std::filesystem::path& target) {
  ACTS_DEBUG("Creating writer for: " << target);

  auto format = HepMC3Util::formatFromFilename(target.string());

  if (format == HepMC3Util::Format::root) {
    ACTS_DEBUG("~> Root");
    if (m_cfg.compression != HepMC3Util::Compression::none) {
      ACTS_ERROR("~~> Compression not supported for Root");
      throw std::invalid_argument("Compression not supported for Root");
    }
#ifdef ACTS_HEPMC3_ROOT_SUPPORT
    return std::make_unique<HepMC3::WriterRootTree>(target);
#else
    ACTS_ERROR("~~> Root support not enabled in HepMC3");
    throw std::runtime_error("Root support not enabled in HepMC3");
#endif
  } else {
    std::filesystem::path path =
        target.string() +
        std::string{HepMC3Util::compressionExtension(m_cfg.compression)};
    ACTS_DEBUG("~> Ascii (=> " << path << ")");

    switch (m_cfg.compression) {
      case HepMC3Util::Compression::none:
        ACTS_DEBUG("~~> uncompressed");
        return std::make_unique<HepMC3::WriterAscii>(path);
#ifdef HEPMC3_USE_COMPRESSION
      case HepMC3Util::Compression::zlib:
        ACTS_DEBUG("~~> GZ");
        return std::make_unique<
            HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::z>>(
            path);
      case HepMC3Util::Compression::lzma:
        ACTS_DEBUG("~~> LZMA");
        return std::make_unique<
            HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::lzma>>(
            path);
      case HepMC3Util::Compression::bzip2:
        ACTS_DEBUG("~~> BZ2");
        return std::make_unique<
            HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::bz2>>(
            path);
      case HepMC3Util::Compression::zstd:
        ACTS_DEBUG("~~> ZSTD");
        return std::make_unique<
            HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::zstd>>(
            path);
#endif
      default:
        throw std::invalid_argument{"Unknown compression value"};
    }
  }
}

HepMC3Writer::~HepMC3Writer() = default;

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

  if (m_cfg.perEvent) {
    std::filesystem::path perEventFile =
        perEventFilepath(m_cfg.outputPath.parent_path(),
                         m_cfg.outputPath.filename().string(), ctx.eventNumber);

    ACTS_VERBOSE("Writing per-event file " << perEventFile);
    auto writer = createWriter(perEventFile);

    writer->write_event(*event);
    auto result = ProcessCode::SUCCESS;
    if (writer->failed()) {
      ACTS_ERROR("Failed to write event number: " << ctx.eventNumber);
      result = ProcessCode::ABORT;
    }
    writer->close();
    return result;
  }

  if (!m_cfg.writeEventsInOrder) {
    std::scoped_lock lock{m_mutex};
    // Unconditionally write events in whatever order they come in
    m_writer->write_event(*event);
    if (m_writer->failed()) {
      ACTS_ERROR("Failed to write event number: " << ctx.eventNumber);
      return ProcessCode::ABORT;
    }
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
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
