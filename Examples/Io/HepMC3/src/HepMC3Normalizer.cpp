// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Normalizer.hpp"

#include <chrono>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderFactory.h>
#include <HepMC3/Writer.h>
#include <HepMC3/WriterAscii.h>

#ifdef HEPMC3_USE_COMPRESSION
#include <HepMC3/CompressedIO.h>
#endif

#ifdef ACTS_HEPMC3_ROOT_SUPPORT
#include <HepMC3/WriterRootTree.h>
#endif

#include <nlohmann/json.hpp>

namespace ActsExamples {

namespace {

/// Write metadata sidecar file
void writeMetadata(const std::filesystem::path& hepmcFile,
                   std::size_t eventCount) {
  auto metaPath = hepmcFile.string() + ".json";

  try {
    nlohmann::json j;
    j["num_events"] = eventCount;

    std::ofstream out(metaPath);
    if (out.is_open()) {
      out << j.dump(2);
    }
  } catch (...) {
    // Silent failure - not critical
  }
}

#ifdef HEPMC3_USE_COMPRESSION
/// Convert HepMC3Util::Compression enum to bxz::Compression
bxz::Compression compressionEnumToBxz(HepMC3Util::Compression compression) {
  using enum HepMC3Util::Compression;
  switch (compression) {
    case none:
      return bxz::Compression::plaintext;
#ifdef HEPMC3_Z_SUPPORT
    case zlib:
      return bxz::Compression::z;
#endif
#ifdef HEPMC3_LZMA_SUPPORT
    case lzma:
      return bxz::Compression::lzma;
#endif
#ifdef HEPMC3_BZ2_SUPPORT
    case bzip2:
      return bxz::Compression::bz2;
#endif
#ifdef HEPMC3_ZSTD_SUPPORT
    case zstd:
      return bxz::Compression::zstd;
#endif
    default:
      throw std::invalid_argument("Unsupported compression type");
  }
}
#endif

/// Create HepMC3 writer with optional compression
std::unique_ptr<HepMC3::Writer> createWriter(
    const std::filesystem::path& path, HepMC3Util::Format format,
    HepMC3Util::Compression compression, int compressionLevel) {
  // ROOT format
  if (format == HepMC3Util::Format::root) {
    if (compression != HepMC3Util::Compression::none) {
      throw std::invalid_argument(
          "ROOT format does not support compression parameter (ROOT has its "
          "own compression)");
    }
#ifdef ACTS_HEPMC3_ROOT_SUPPORT
    return std::make_unique<HepMC3::WriterRootTree>(path.string());
#else
    throw std::runtime_error("ROOT support not enabled in this build");
#endif
  }

  // ASCII format
#ifdef HEPMC3_USE_COMPRESSION
  bxz::Compression comp = compressionEnumToBxz(compression);

  if (comp == bxz::Compression::plaintext) {
    // Uncompressed ASCII
    return std::make_unique<HepMC3::WriterAscii>(path.string());
  } else {
    // Compressed ASCII with custom level
    auto stream =
        std::make_shared<bxz::ofstream>(path.string(), comp, compressionLevel);
    return std::make_unique<HepMC3::WriterAscii>(stream);
  }
#else
  if (compression != HepMC3Util::Compression::none) {
    throw std::runtime_error("Compression support not enabled in this build");
  }
  return std::make_unique<HepMC3::WriterAscii>(path.string());
#endif
}

/// Generate output filename based on parameters
std::string generateOutputFilename(const std::string& prefix,
                                   HepMC3Util::Format format,
                                   HepMC3Util::Compression compression,
                                   std::size_t eventIndex) {
  // Multi-file mode: append event index
  std::string filename = std::format("{}_{:06d}", prefix, eventIndex);

  // Add format extension
  if (format == HepMC3Util::Format::ascii) {
    filename += ".hepmc3";
    filename += std::string(HepMC3Util::compressionExtension(compression));
  } else {
    filename += ".root";
  }

  return filename;
}

/// Format file size as human-readable string
std::string formatSize(std::uintmax_t bytes) {
  const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
  int unitIndex = 0;
  double size = static_cast<double>(bytes);

  while (size >= 1024.0 && unitIndex < 4) {
    size /= 1024.0;
    unitIndex++;
  }

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << size << " " << units[unitIndex];
  return oss.str();
}

}  // namespace

HepMC3Normalizer::HepMC3Normalizer(Config cfg) : m_cfg(std::move(cfg)) {
  // Validate configuration
  if (m_cfg.inputFiles.empty()) {
    throw std::invalid_argument("No input files specified");
  }

  if (!m_cfg.singleOutputPath && m_cfg.eventsPerFile == 0) {
    throw std::invalid_argument(
        "events-per-file must be > 0 in multi-file mode");
  }

  if (m_cfg.compressionLevel < 0 || m_cfg.compressionLevel > 19) {
    throw std::invalid_argument("compression-level must be 0-19");
  }

  if (m_cfg.format == HepMC3Util::Format::root &&
      m_cfg.compression != HepMC3Util::Compression::none) {
    throw std::invalid_argument(
        "ROOT format does not support compression parameter");
  }
}

HepMC3Normalizer::Result HepMC3Normalizer::normalize() {
  Result result;

  // Determine output directory
  std::filesystem::path outputDir;
  if (m_cfg.singleOutputPath) {
    if (m_cfg.singleOutputPath->has_parent_path()) {
      outputDir = m_cfg.singleOutputPath->parent_path();
    } else {
      outputDir = ".";
    }
  } else {
    outputDir = m_cfg.outputDir;
  }

  // Create output directory
  std::filesystem::create_directories(outputDir);

  if (m_cfg.verbose) {
    std::cerr << "Configuration:\n";
    std::cerr << "  Input files: " << m_cfg.inputFiles.size() << "\n";
    if (m_cfg.singleOutputPath) {
      std::cerr << "  Single output file: " << *m_cfg.singleOutputPath << "\n";
      std::cerr << "  Format: " << m_cfg.format << "\n";
      std::cerr << "  Compression: " << m_cfg.compression << " (level "
                << m_cfg.compressionLevel << ")\n";
    } else {
      std::cerr << "  Output dir: " << outputDir << "\n";
      std::cerr << "  Output prefix: " << m_cfg.outputPrefix << "\n";
      std::cerr << "  Events per file: " << m_cfg.eventsPerFile << "\n";
      std::cerr << "  Format: " << m_cfg.format << "\n";
      std::cerr << "  Compression: " << m_cfg.compression << " (level "
                << m_cfg.compressionLevel << ")\n";
    }
    if (m_cfg.maxEvents > 0) {
      std::cerr << "  Max events to read: " << m_cfg.maxEvents << "\n";
    }
    std::cerr << "\n";
  }

  // Processing state
  std::size_t globalEventIndex = 0;
  std::size_t eventsInCurrentFile = 0;
  std::unique_ptr<HepMC3::Writer> currentWriter;
  std::filesystem::path currentOutputPath;

  // Process each input file
  for (const auto& inputFile : m_cfg.inputFiles) {
    if (m_cfg.verbose) {
      std::cerr << "Reading " << inputFile << "...\n";
    }

    if (!std::filesystem::exists(inputFile)) {
      std::cerr << "WARNING: File not found: " << inputFile << "\n";
      continue;
    }

    // Track input file size
    result.totalInputSize += std::filesystem::file_size(inputFile);

    auto reader = HepMC3::deduce_reader(inputFile.string());
    if (!reader) {
      std::cerr << "ERROR: Failed to open " << inputFile << "\n";
      continue;
    }

    HepMC3::GenEvent event;
    std::size_t eventsReadFromFile = 0;
    std::size_t lastProgressUpdate = 0;
    constexpr std::size_t progressInterval = 100;

    while (!reader->failed()) {
      // Check if we've reached the maximum number of events
      if (m_cfg.maxEvents > 0 && globalEventIndex >= m_cfg.maxEvents) {
        break;
      }

      auto readStart = std::chrono::high_resolution_clock::now();
      reader->read_event(event);
      auto readEnd = std::chrono::high_resolution_clock::now();
      result.totalReadTime +=
          std::chrono::duration<double>(readEnd - readStart).count();

      if (reader->failed()) {
        break;
      }

      // Show progress
      if (m_cfg.verbose &&
          (eventsReadFromFile - lastProgressUpdate) >= progressInterval) {
        std::cerr << "    Progress: " << eventsReadFromFile
                  << " events read...\r" << std::flush;
        lastProgressUpdate = eventsReadFromFile;
      }

      // Create new output file if needed
      if (eventsInCurrentFile == 0) {
        // Close previous file
        if (currentWriter) {
          currentWriter->close();

          // Get file size and write metadata
          if (std::filesystem::exists(currentOutputPath)) {
            auto fileSize = std::filesystem::file_size(currentOutputPath);
            result.totalOutputSize += fileSize;

            writeMetadata(currentOutputPath,
                          m_cfg.singleOutputPath ? globalEventIndex
                                                 : m_cfg.eventsPerFile);
            result.outputFiles.push_back(currentOutputPath);

            if (m_cfg.verbose) {
              std::size_t events =
                  m_cfg.singleOutputPath ? globalEventIndex : m_cfg.eventsPerFile;
              double sizePerEvent =
                  static_cast<double>(fileSize) / static_cast<double>(events);

              std::cerr << "  Wrote " << currentOutputPath << " (" << events
                        << " events, " << formatSize(fileSize) << ", "
                        << formatSize(static_cast<std::uintmax_t>(sizePerEvent))
                        << "/event)\n";
            }
          }
        }

        // Generate output path based on mode
        if (m_cfg.singleOutputPath) {
          currentOutputPath = *m_cfg.singleOutputPath;
        } else {
          std::string filename = generateOutputFilename(
              m_cfg.outputPrefix, m_cfg.format, m_cfg.compression,
              globalEventIndex);
          currentOutputPath = outputDir / filename;
        }
        currentWriter = createWriter(currentOutputPath, m_cfg.format,
                                     m_cfg.compression, m_cfg.compressionLevel);
      }

      // Set event number
      event.set_event_number(globalEventIndex);

      // Clear run info from events that are not the first in their output file
      if (eventsInCurrentFile > 0) {
        event.set_run_info(nullptr);
      }

      auto writeStart = std::chrono::high_resolution_clock::now();
      currentWriter->write_event(event);
      auto writeEnd = std::chrono::high_resolution_clock::now();
      result.totalWriteTime +=
          std::chrono::duration<double>(writeEnd - writeStart).count();

      globalEventIndex++;
      eventsInCurrentFile++;
      eventsReadFromFile++;

      // Close file if chunk is complete (only in multi-file mode)
      if (!m_cfg.singleOutputPath &&
          eventsInCurrentFile >= m_cfg.eventsPerFile) {
        currentWriter->close();

        // Get file size
        if (std::filesystem::exists(currentOutputPath)) {
          auto fileSize = std::filesystem::file_size(currentOutputPath);
          result.totalOutputSize += fileSize;

          writeMetadata(currentOutputPath, eventsInCurrentFile);
          result.outputFiles.push_back(currentOutputPath);

          if (m_cfg.verbose) {
            double sizePerEvent = static_cast<double>(fileSize) /
                                  static_cast<double>(eventsInCurrentFile);

            std::cerr << "  Wrote " << currentOutputPath << " ("
                      << eventsInCurrentFile << " events, "
                      << formatSize(fileSize) << ", "
                      << formatSize(static_cast<std::uintmax_t>(sizePerEvent))
                      << "/event)\n";
          }
        }

        currentWriter.reset();
        eventsInCurrentFile = 0;
      }
    }

    reader->close();

    if (m_cfg.verbose) {
      // Clear progress line and print final count
      std::cerr << "  Read " << eventsReadFromFile << " events from "
                << inputFile << "                    \n";
    }

    // Check if we've reached the maximum number of events
    if (m_cfg.maxEvents > 0 && globalEventIndex >= m_cfg.maxEvents) {
      if (m_cfg.verbose) {
        std::cerr << "Reached maximum event limit (" << m_cfg.maxEvents
                  << ")\n";
      }
      break;
    }
  }

  // Close final file
  if (currentWriter && eventsInCurrentFile > 0) {
    currentWriter->close();

    // Get file size
    if (std::filesystem::exists(currentOutputPath)) {
      auto fileSize = std::filesystem::file_size(currentOutputPath);
      result.totalOutputSize += fileSize;

      writeMetadata(currentOutputPath, eventsInCurrentFile);
      result.outputFiles.push_back(currentOutputPath);

      if (m_cfg.verbose) {
        double sizePerEvent = static_cast<double>(fileSize) /
                              static_cast<double>(eventsInCurrentFile);

        std::cerr << "  Wrote " << currentOutputPath << " ("
                  << eventsInCurrentFile << " events, "
                  << formatSize(fileSize) << ", "
                  << formatSize(static_cast<std::uintmax_t>(sizePerEvent))
                  << "/event)\n";
      }
    }
  }

  result.numEvents = globalEventIndex;

  // Print summary
  if (m_cfg.verbose) {
    std::cerr << "\nSummary:\n";
    if (m_cfg.singleOutputPath) {
      std::cerr << "  Processed " << globalEventIndex
                << " events into single file\n";
    } else {
      std::size_t totalFiles =
          (globalEventIndex + m_cfg.eventsPerFile - 1) / m_cfg.eventsPerFile;
      std::cerr << "  Processed " << globalEventIndex << " events into "
                << totalFiles << " file(s)\n";
    }

    if (globalEventIndex > 0) {
      double bytesPerEvent = static_cast<double>(result.totalOutputSize) /
                             static_cast<double>(globalEventIndex);
      std::cerr << "  Total input size:  " << formatSize(result.totalInputSize)
                << "\n";
      std::cerr << "  Total output size: " << formatSize(result.totalOutputSize)
                << " (" << formatSize(static_cast<std::uintmax_t>(bytesPerEvent))
                << "/event)\n";

      if (result.totalInputSize > 0) {
        double ratio = static_cast<double>(result.totalOutputSize) /
                       static_cast<double>(result.totalInputSize);
        std::cerr << "  Compression ratio: " << std::fixed
                  << std::setprecision(2) << (ratio * 100.0) << "%\n";
      }
    } else {
      std::cerr << "  Total input size:  " << formatSize(result.totalInputSize)
                << "\n";
      std::cerr << "  Total output size: " << formatSize(result.totalOutputSize)
                << "\n";
    }

    // Print timing information
    if (globalEventIndex > 0) {
      std::cerr << "\nTiming breakdown:\n";
      std::cerr << "  Reading events:  " << std::fixed << std::setprecision(3)
                << result.totalReadTime << " s ("
                << (result.totalReadTime / globalEventIndex * 1000.0)
                << " ms/event)\n";
      std::cerr << "  Writing events:  " << std::fixed << std::setprecision(3)
                << result.totalWriteTime << " s ("
                << (result.totalWriteTime / globalEventIndex * 1000.0)
                << " ms/event)\n";

      double totalProcessingTime = result.totalReadTime + result.totalWriteTime;
      std::cerr << "  Total processing: " << std::fixed << std::setprecision(3)
                << totalProcessingTime << " s\n";

      // Show percentage breakdown
      if (totalProcessingTime > 0) {
        std::cerr << "\nTime distribution:\n";
        std::cerr << "  Reading: " << std::fixed << std::setprecision(1)
                  << (result.totalReadTime / totalProcessingTime * 100.0)
                  << "%\n";
        std::cerr << "  Writing: " << std::fixed << std::setprecision(1)
                  << (result.totalWriteTime / totalProcessingTime * 100.0)
                  << "%\n";
      }
    }
  }

  return result;
}

}  // namespace ActsExamples
