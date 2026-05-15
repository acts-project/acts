// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

#include <filesystem>
#include <map>
#include <mutex>
#include <semaphore>
#include <string>

namespace HepMC3 {
class GenEvent;
class Writer;
}  // namespace HepMC3

namespace ActsExamples {

/// HepMC3 event writer.
class HepMC3Writer final : public WriterT<std::shared_ptr<HepMC3::GenEvent>> {
 public:
  struct Config {
    /// The output file path for writing HepMC3 events.
    /// All events will be written to a single file.
    std::filesystem::path outputPath;

    /// The input collection
    std::string inputEvent;

    /// The compression mode to use for the output file.
    /// If set, the compression is explicitly specified and the outputPath can
    /// omit the compression extension (e.g., ".hepmc3" instead of
    /// ".hepmc3.gz"). If both compression and the extension in outputPath are
    /// specified, they must be consistent.
    /// If not set (nullopt), the compression is deduced from the outputPath
    /// extension.
    std::optional<HepMC3Util::Compression> compression = std::nullopt;

    /// Write events in order of their event number. If `false`, events are
    /// written in whatever order they are processed.
    bool writeEventsInOrder = true;

    /// The maximum number of events to keep in the queue before writing
    /// The value depends on the number of core and therefore the worst-case
    /// distance between the *current* and the next event we can write to write
    /// in order.
    /// @note If @c writeEventsInOrder is false, this value is ignored.
    unsigned long maxEventsPending = 128;
  };

  /// Construct the writer.
  ///
  /// @param [in] config Config of the writer
  /// @param [in] level The level of the logger
  HepMC3Writer(const Config& config, Acts::Logging::Level level);

  ~HepMC3Writer() override;

  ProcessCode beginEvent(std::size_t threadId) override;

  /// Writing events to file.
  ///
  /// @param [in] ctx The context of this algorithm
  /// @param [in] event The recorded HepMC3 event
  ///
  /// @return Code describing whether the writing was successful
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::shared_ptr<HepMC3::GenEvent>& event) override;

  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Resolve and validate compression settings.
  /// Returns the resolved compression mode and the base path (without
  /// compression extension).
  std::pair<HepMC3Util::Compression, std::filesystem::path> resolveCompression(
      const std::filesystem::path& path) const;

  /// The configuration of this writer
  Config m_cfg;

  /// The resolved compression mode (deduced from config or path)
  HepMC3Util::Compression m_compression;

  /// The actual output path with compression extension
  std::filesystem::path m_actualOutputPath;

  std::mutex m_mutex;

  std::size_t m_nextEvent = 0;
  std::vector<std::pair<long long, std::shared_ptr<HepMC3::GenEvent>>>
      m_eventQueue;

  // For reporting at the end of the job
  std::size_t m_maxEventQueueSize = 0;

  std::counting_semaphore<> m_queueSemaphore;
  std::atomic<unsigned int> m_waiting = 0;

  std::unique_ptr<HepMC3::Writer> m_writer;

  // Track the number of events written for metadata
  std::size_t m_eventsWritten = 0;
};

}  // namespace ActsExamples
