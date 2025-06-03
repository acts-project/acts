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
    /// If true, one file per event is written with the event number appended to
    /// the filename
    bool perEvent = false;

    /// The output file path for writing HepMC3 events.
    ///
    /// This path is handled differently based on the perEvent flag:
    /// - If perEvent is false: The path points to a single file where all
    /// events will be written
    /// - If perEvent is true: The path is used as a template for creating
    /// per-event files
    ///   in the format "event{number}-{filename}" in the parent directory
    ///
    /// When in per-event mode, the writer uses perEventFilepath() to generate
    /// the appropriate filename for each event.
    std::filesystem::path outputPath;

    /// The input collection
    std::string inputEvent;

    /// The compression mode to use for the output file
    HepMC3Util::Compression compression = HepMC3Util::Compression::none;

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
  ProcessCode writeT(const ActsExamples::AlgorithmContext& ctx,
                     const std::shared_ptr<HepMC3::GenEvent>& event) override;

  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  std::unique_ptr<HepMC3::Writer> createWriter(
      const std::filesystem::path& target);

  /// The configuration of this writer
  Config m_cfg;

  std::mutex m_mutex;

  std::size_t m_nextEvent = 0;
  std::vector<std::pair<long long, std::shared_ptr<HepMC3::GenEvent>>>
      m_eventQueue;

  // For reporting at the end of the job
  std::size_t m_maxEventQueueSize = 0;

  std::counting_semaphore<> m_queueSemaphore;
  std::atomic<unsigned int> m_waiting = 0;

  std::unique_ptr<HepMC3::Writer> m_writer;
};

}  // namespace ActsExamples
