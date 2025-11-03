// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <vector>

namespace HepMC3 {
class GenEvent;
class Reader;
class Writer;
}  // namespace HepMC3

namespace Acts {
class Logger;
}

namespace ActsExamples::HepMC3Util {

void mergeEvents(HepMC3::GenEvent& event,
                 std::span<const HepMC3::GenEvent*> genEvents,
                 const Acts::Logger& logger);

void mergeEvents(HepMC3::GenEvent& event,
                 std::span<std::shared_ptr<const HepMC3::GenEvent>> genEvents,
                 const Acts::Logger& logger);

enum class Compression { none, zlib, lzma, bzip2, zstd };

std::ostream& operator<<(std::ostream& os, HepMC3Util::Compression compression);

std::span<const Compression> availableCompressionModes();

std::string_view compressionExtension(Compression compression);

enum class Format { ascii, root };

std::ostream& operator<<(std::ostream& os, Format format);

std::span<const Format> availableFormats();

Format formatFromFilename(std::string_view filename);

/// Result of HepMC3 file normalization
struct NormalizeResult {
  /// Number of events processed
  std::size_t numEvents = 0;

  /// Output files created
  std::vector<std::filesystem::path> outputFiles;

  /// Total input size in bytes
  std::size_t totalInputSize = 0;

  /// Total output size in bytes
  std::size_t totalOutputSize = 0;

  /// Total time spent reading (seconds)
  double totalReadTime = 0.0;

  /// Total time spent writing (seconds)
  double totalWriteTime = 0.0;
};

/// Normalize and optionally chunk HepMC3 files.
///
/// Reads one or more HepMC3 files, normalizes event numbers, and writes them
/// to output files. Can write to a single output file or chunk events into
/// multiple files.
///
/// @param inputFiles Input HepMC files to normalize
/// @param singleOutputPath Single output file path (optional). If specified,
///        all events are written to this file. Format and compression are
///        auto-detected from filename. Mutually exclusive with chunking.
/// @param outputDir Output directory for multi-file mode
/// @param outputPrefix Output file prefix for multi-file mode
/// @param eventsPerFile Number of events per output file in multi-file mode
/// @param maxEvents Maximum number of events to process (`std::nullopt` = all events)
/// @param format Output format (ascii or root)
/// @param compression Compression type
/// @param compressionLevel Compression level (0-19, higher = more compression)
/// @param verbose Enable verbose output
/// @return Result with statistics and list of created files
NormalizeResult normalizeFiles(
    const std::vector<std::filesystem::path>& inputFiles,
    std::optional<std::filesystem::path> singleOutputPath = std::nullopt,
    const std::filesystem::path& outputDir = ".",
    const std::string& outputPrefix = "events",
    std::optional<std::size_t> eventsPerFile = 10000,
    std::optional<std::size_t> maxEvents = std::nullopt,
    Format format = Format::ascii, Compression compression = Compression::none,
    int compressionLevel = 6, bool verbose = false);

/// Wrapper around HepMC3::deduce_reader to isolate problematic HepMC3 headers
/// that have multiple definition issues in versions < 3.3.0
///
/// @param filename Path to the HepMC3 file
/// @return Shared pointer to the appropriate HepMC3 reader
std::shared_ptr<HepMC3::Reader> deduceReader(const std::string& filename);

/// Wrapper around HepMC3 writer creation to isolate problematic HepMC3 headers
/// that have multiple definition issues in versions < 3.3.0
///
/// @param path Path to the output file
/// @param format Output format (ascii or root)
/// @param compression Compression type
/// @return Unique pointer to the appropriate HepMC3 writer
std::unique_ptr<HepMC3::Writer> createWriter(const std::filesystem::path& path,
                                             Format format,
                                             Compression compression);

}  // namespace ActsExamples::HepMC3Util
