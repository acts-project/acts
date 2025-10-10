// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace ActsExamples {

/// HepMC3 file normalizer and chunker.
/// Reads HepMC3 files, normalizes event numbers, and optionally chunks
/// into multiple output files.
class HepMC3Normalizer {
 public:
  /// Configuration for the normalizer
  struct Config {
    /// Input HepMC files to normalize
    std::vector<std::filesystem::path> inputFiles;

    /// Single output file path (optional - mutually exclusive with chunking)
    std::optional<std::filesystem::path> singleOutputPath;

    /// Output directory for multi-file mode
    std::filesystem::path outputDir = ".";

    /// Output file prefix for multi-file mode
    std::string outputPrefix = "events";

    /// Number of events per output file (0 = single file mode)
    std::size_t eventsPerFile = 10000;

    /// Maximum number of events to process (0 = all)
    std::size_t maxEvents = 0;

    /// Output format (ascii or root)
    HepMC3Util::Format format = HepMC3Util::Format::ascii;

    /// Compression type
    HepMC3Util::Compression compression = HepMC3Util::Compression::none;

    /// Compression level (0-19)
    int compressionLevel = 6;

    /// Enable verbose output
    bool verbose = false;
  };

  /// Result of normalization
  struct Result {
    /// Number of events processed
    std::size_t numEvents = 0;

    /// Output files created
    std::vector<std::filesystem::path> outputFiles;

    /// Total input size in bytes
    std::uintmax_t totalInputSize = 0;

    /// Total output size in bytes
    std::uintmax_t totalOutputSize = 0;

    /// Total time spent reading (seconds)
    double totalReadTime = 0.0;

    /// Total time spent writing (seconds)
    double totalWriteTime = 0.0;
  };

  /// Construct normalizer with configuration
  explicit HepMC3Normalizer(Config cfg);

  /// Run the normalization
  Result normalize();

  /// Get the configuration
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
