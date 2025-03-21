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
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <filesystem>
#include <string>

namespace HepMC3 {
class GenEvent;
class Writer;
}  // namespace HepMC3

namespace ActsExamples {

/// HepMC3 event writer.
class HepMC3AsciiWriter final : public WriterT<std::vector<HepMC3::GenEvent>> {
 public:
  struct Config {
    /// If true, one file per event is written with the event number appended to
    /// the filename
    bool perEvent = false;

    /// The output file path
    std::filesystem::path outputPath;

    // The input collection
    std::string inputEvents;
  };

  /// Construct the writer.
  ///
  /// @param [in] config Config of the writer
  /// @param [in] level The level of the logger
  HepMC3AsciiWriter(const Config& config, Acts::Logging::Level level);

  ~HepMC3AsciiWriter();

  /// Writing events to file.
  ///
  /// @param [in] ctx The context of this algorithm
  /// @param [in] events The recorded HepMC3 events
  ///
  /// @return Code describing whether the writing was successful
  ProcessCode writeT(const ActsExamples::AlgorithmContext& ctx,
                     const std::vector<HepMC3::GenEvent>& events) override;

  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// The configuration of this writer
  Config m_cfg;

  std::mutex m_mutex;

  std::unique_ptr<HepMC3::Writer> m_writer;
};

}  // namespace ActsExamples
