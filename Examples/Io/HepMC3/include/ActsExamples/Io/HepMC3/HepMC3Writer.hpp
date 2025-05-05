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
#include <mutex>
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
  };

  /// Construct the writer.
  ///
  /// @param [in] config Config of the writer
  /// @param [in] level The level of the logger
  HepMC3Writer(const Config& config, Acts::Logging::Level level);

  ~HepMC3Writer() override;

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

  std::unique_ptr<HepMC3::Writer> m_writer;
};

}  // namespace ActsExamples
