// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

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
class GenParticle;
class GenVertex;
class FourVector;
class Reader;
class Writer;
}  // namespace HepMC3

namespace Acts {
class Logger;
}

namespace ActsExamples::HepMC3Util {

/// Merge multiple generator events into a single target event.
///
/// @param event Destination event that receives merged particles
/// @param genEvents Collection of events to merge into @p event
/// @param logger Logger used for diagnostic messages
void mergeEvents(HepMC3::GenEvent& event,
                 std::span<const HepMC3::GenEvent*> genEvents,
                 const Acts::Logger& logger);

/// Merge multiple generator events into a single target event.
///
/// @param event Destination event that receives merged particles
/// @param genEvents Collection of shared events to merge into @p event
/// @param logger Logger used for diagnostic messages
void mergeEvents(HepMC3::GenEvent& event,
                 std::span<std::shared_ptr<const HepMC3::GenEvent>> genEvents,
                 const Acts::Logger& logger);

constexpr int kBeamParticleStatus = 4;
constexpr int kUndecayedParticleStatus = 1;
constexpr int kDecayedParticleStatus = 2;

/// Supported compression codecs for HepMC3 files.
enum class Compression { none, zlib, lzma, bzip2, zstd };

/// Stream operator for human-readable `Compression` output.
///
/// @param os Output stream receiving the textual representation
/// @param compression Compression mode to serialize
/// @return Reference to @p os for chaining
std::ostream& operator<<(std::ostream& os, HepMC3Util::Compression compression);

/// List of compression modes available in the linked HepMC3 build.
///
/// @return Span of supported compression values
std::span<const Compression> availableCompressionModes();

/// File extension associated with a compression mode (e.g. `.gz`).
///
/// @param compression Compression mode to query
/// @return Extension string including leading dot when present, empty otherwise
std::string_view compressionExtension(Compression compression);

/// Deduce compression mode from a file path.
///
/// @param filename Path used to infer compression from its suffix
/// @return Compression inferred from @p filename
Compression compressionFromFilename(const std::filesystem::path& filename);

/// Supported HepMC3 output formats.
enum class Format { ascii, root };

/// Stream operator for human-readable `Format` output.
///
/// @param os Output stream receiving the textual representation
/// @param format File format to serialize
/// @return Reference to @p os for chaining
std::ostream& operator<<(std::ostream& os, Format format);

/// List of output formats available in the linked HepMC3 build.
///
/// @return Span of supported format values
std::span<const Format> availableFormats();

/// Deduce HepMC3 format from a file path.
///
/// @param filename Path used to infer format from its suffix
/// @return Format inferred from @p filename
Format formatFromFilename(const std::filesystem::path& filename);

static constexpr std::string_view kEventGeneratorIndexAttribute =
    "acts_gen_event_index";

int eventGeneratorIndex(const HepMC3::GenParticle& particle);
int eventGeneratorIndex(const HepMC3::GenVertex& vertex);

Acts::Vector4 convertPosition(const HepMC3::FourVector& vec);

std::vector<const HepMC3::GenVertex*> findHardScatterVertices(
    const HepMC3::GenEvent& event);

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
