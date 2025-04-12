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

#include <boost/algorithm/string/join.hpp>

namespace ActsExamples {

HepMC3Writer::HepMC3Writer(const Config& config, Acts::Logging::Level level)
    : WriterT(config.inputEvent, "HepMC3Writer", level), m_cfg(config) {
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
}

std::unique_ptr<HepMC3::Writer> HepMC3Writer::createWriter(
    const std::filesystem::path& target) {
  std::filesystem::path path =
      target.string() +
      std::string{HepMC3Util::compressionExtension(m_cfg.compression)};

  switch (m_cfg.compression) {
    case HepMC3Util::Compression::none:
      return std::make_unique<HepMC3::WriterAscii>(path);
#ifdef HEPMC3_USE_COMPRESSION
    case HepMC3Util::Compression::zlib:
      return std::make_unique<
          HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::z>>(path);
    case HepMC3Util::Compression::lzma:
      return std::make_unique<
          HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::lzma>>(
          path);
    case HepMC3Util::Compression::bzip2:
      return std::make_unique<
          HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::bz2>>(
          path);
    case HepMC3Util::Compression::zstd:
      return std::make_unique<
          HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::zstd>>(
          path);
#endif
    default:
      throw std::invalid_argument{"Unknown compression value"};
  }
}

HepMC3Writer::~HepMC3Writer() = default;

ProcessCode HepMC3Writer::writeT(
    const AlgorithmContext& ctx,
    const std::shared_ptr<HepMC3::GenEvent>& event) {
  ACTS_VERBOSE("Writing " << event->particles().size() << " particles to "
                          << m_cfg.outputPath);

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

  // Lock is needed both for the queueing as well as the flushing
  std::scoped_lock lock(m_mutex);

  if (auto pc = queueForWriting(ctx.eventNumber, event);
      pc != ProcessCode::SUCCESS) {
    return pc;
  }

  flushQueue();

  return ProcessCode::SUCCESS;
}

ProcessCode HepMC3Writer::queueForWriting(
    std::size_t eventNumber, std::shared_ptr<HepMC3::GenEvent> event) {
  ACTS_DEBUG("Queueing event_number=" << eventNumber << ", current_length="
                                      << m_eventQueue.size());
  if (m_eventQueue.size() > m_cfg.maxEventsPending) {
    ACTS_ERROR("Queue size of "
               << m_eventQueue.size() << " would exceed maximum of "
               << m_cfg.maxEventsPending
               << ". Cannot proceed without changing event ordering");
    return ProcessCode::ABORT;
  }

  ACTS_VERBOSE("Finding insert location for event number: " << eventNumber);
  auto it = std::ranges::upper_bound(m_eventQueue, eventNumber, {},
                                     [](const auto& v) { return v.first; });
  if (it == m_eventQueue.end()) {
    ACTS_VERBOSE("Insert location is at the end of the queue");
  } else {
    ACTS_VERBOSE("Insert location is before event number: " << it->first);
  }

  m_eventQueue.insert(it, {eventNumber, std::move(event)});

  m_maxEventQueueSize = std::max(m_maxEventQueueSize, m_eventQueue.size());

  return ProcessCode::SUCCESS;
}

ProcessCode HepMC3Writer::flushQueue() {
  ACTS_DEBUG("Flushing queue, next_event=" << m_nextEvent << ", queue_length="
                                           << m_eventQueue.size());

  ACTS_VERBOSE("queue=[" << [&]() {
    std::vector<std::string> numbers;
    numbers.reserve(m_eventQueue.size());
    std::ranges::transform(
        m_eventQueue, std::back_inserter(numbers),
        [](const auto& pair) { return std::to_string(pair.first); });

    return boost::algorithm::join(numbers, ", ");
  }() << "]");

  while (!m_eventQueue.empty() && m_eventQueue.front().first == m_nextEvent) {
    auto next = std::move(m_eventQueue.front());
    ACTS_VERBOSE("Writing event number: " << next.first);
    m_eventQueue.erase(m_eventQueue.begin());

    m_writer->write_event(*next.second);
    if (m_writer->failed()) {
      ACTS_ERROR("Failed to write event number: " << next.first);
      return ProcessCode::ABORT;
    }
    m_nextEvent++;
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
