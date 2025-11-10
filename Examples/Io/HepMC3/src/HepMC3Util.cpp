// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Reader.h>
#include <HepMC3/Version.h>
#include <HepMC3/Writer.h>
#include <nlohmann/json.hpp>

namespace ActsExamples {

namespace {
template <typename T>
void mergeEventsImpl(HepMC3::GenEvent& event, std::span<T> genEvents,
                     const Acts::Logger& logger) {
  Acts::AveragingScopedTimer mergeTimer("Merging HepMC3 events", logger(),
                                        Acts::Logging::DEBUG);

  std::vector<std::shared_ptr<HepMC3::GenParticle>> particles;

  // Loop once to find the total size we'll need
  std::size_t nParticles = 0;
  std::size_t nVertices = 0;
  for (const auto& genEvent : genEvents) {
    nParticles += genEvent->particles().size();
    nVertices += genEvent->vertices().size();
  }

  event.reserve(nParticles, nVertices);

  for (const auto& genEvent : genEvents) {
    auto sample = mergeTimer.sample();
    particles.clear();
    particles.reserve(genEvent->particles().size());

    auto copyAttributes = [&](const auto& src, auto& dst) {
      for (const auto& attr : src.attribute_names()) {
        auto value = src.attribute_as_string(attr);
        dst.add_attribute(attr,
                          std::make_shared<HepMC3::StringAttribute>(value));
      }
    };

    copyAttributes(*genEvent, event);

    // Add to combined event
    for (const auto& srcParticle : genEvent->particles()) {
      if (srcParticle->id() - 1 != static_cast<int>(particles.size())) {
        throw std::runtime_error("Particle id is not consecutive");
      }
      auto particle = std::make_shared<HepMC3::GenParticle>();
      particle->set_momentum(srcParticle->momentum());
      particle->set_generated_mass(srcParticle->generated_mass());
      particle->set_pid(srcParticle->pid());
      particle->set_status(srcParticle->status());

      particles.push_back(particle);
      event.add_particle(particle);

      copyAttributes(*srcParticle, *particle);
    }

    for (const auto& srcVertex : genEvent->vertices()) {
      auto vertex = std::make_shared<HepMC3::GenVertex>(srcVertex->position());
      vertex->set_status(srcVertex->status());

      event.add_vertex(vertex);

      copyAttributes(*srcVertex, *vertex);

      for (const auto& srcParticle : srcVertex->particles_in()) {
        const auto& particle = particles.at(srcParticle->id() - 1);
        vertex->add_particle_in(particle);
      }
      for (const auto& srcParticle : srcVertex->particles_out()) {
        const auto& particle = particles.at(srcParticle->id() - 1);
        vertex->add_particle_out(particle);
      }
    }
  }
}
}  // namespace

void HepMC3Util::mergeEvents(HepMC3::GenEvent& event,
                             std::span<const HepMC3::GenEvent*> genEvents,
                             const Acts::Logger& logger) {
  mergeEventsImpl(event, genEvents, logger);
}

void HepMC3Util::mergeEvents(
    HepMC3::GenEvent& event,
    std::span<std::shared_ptr<const HepMC3::GenEvent>> genEvents,
    const Acts::Logger& logger) {
  mergeEventsImpl(event, genEvents, logger);
}

std::string_view HepMC3Util::compressionExtension(Compression compression) {
  switch (compression) {
    using enum Compression;
    case none:
      return "";
    case zlib:
      return ".gz";
    case lzma:
      return ".xz";
    case bzip2:
      return ".bz2";
    case zstd:
      return ".zst";
    default:
      throw std::invalid_argument{"Unknown compression value"};
  }
}

HepMC3Util::Compression HepMC3Util::compressionFromFilename(
    const std::filesystem::path& filename) {
  using enum Compression;

  std::string filenameStr = filename.string();

  // Check for compression extensions in order
  if (filenameStr.ends_with(".gz")) {
    return zlib;
  }
  if (filenameStr.ends_with(".xz")) {
    return lzma;
  }
  if (filenameStr.ends_with(".bz2")) {
    return bzip2;
  }
  if (filenameStr.ends_with(".zst")) {
    return zstd;
  }

  // No compression extension found
  return none;
}

std::span<const HepMC3Util::Compression>
HepMC3Util::availableCompressionModes() {
  using enum Compression;
  static const auto values = []() -> std::vector<HepMC3Util::Compression> {
    return {
        none,
#ifdef HEPMC3_Z_SUPPORT
        zlib,
#endif
#ifdef HEPMC3_LZMA_SUPPORT
        lzma,
#endif
#ifdef HEPMC3_BZ2_SUPPORT
        bzip2,
#endif
#ifdef HEPMC3_ZSTD_SUPPORT
        zstd,
#endif
    };
  }();
  return values;
}

std::ostream& HepMC3Util::operator<<(std::ostream& os,
                                     HepMC3Util::Compression compression) {
  switch (compression) {
    using enum HepMC3Util::Compression;
    case none:
      return os << "none";
    case zlib:
      return os << "zlib";
    case lzma:
      return os << "lzma";
    case bzip2:
      return os << "bzip2";
    case zstd:
      return os << "zstd";
    default:
      throw std::invalid_argument{"Unknown compression value"};
  }
}

std::ostream& HepMC3Util::operator<<(std::ostream& os,
                                     HepMC3Util::Format format) {
  switch (format) {
    using enum HepMC3Util::Format;
    case ascii:
      return os << "ascii";
    case root:
      return os << "root";
    default:
      throw std::invalid_argument{"Unknown format value"};
  }
}

std::span<const HepMC3Util::Format> HepMC3Util::availableFormats() {
  using enum Format;
  static const auto values = []() -> std::vector<HepMC3Util::Format> {
    return {
        ascii,
#ifdef HEPMC3_ROOT_SUPPORT
        root,
#endif
    };
  }();
  return values;
}

HepMC3Util::Format HepMC3Util::formatFromFilename(
    const std::filesystem::path& filename) {
  using enum Format;

  for (auto compression : availableCompressionModes()) {
    auto ext = compressionExtension(compression);

    if (filename.string().ends_with(".hepmc3" + std::string(ext)) ||
        filename.string().ends_with(".hepmc" + std::string(ext))) {
      return ascii;
    }
  }
  if (filename.string().ends_with(".root")) {
    return root;
  }

  throw std::invalid_argument{"Unknown format extension: " +
                              std::string{filename}};
}

namespace {

// Helper function to format file sizes
std::string formatSize(std::uintmax_t bytes) {
  const std::array<const char*, 5> units = {"B", "KB", "MB", "GB", "TB"};
  std::size_t unit = 0;
  double size = static_cast<double>(bytes);

  while (size >= 1024.0 && unit < units.size() - 1) {
    size /= 1024.0;
    unit++;
  }

  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << size << " " << units[unit];
  return oss.str();
}

// Helper function to write metadata sidecar file
void writeMetadata(const std::filesystem::path& dataFile,
                   std::size_t numEvents) {
  std::filesystem::path metadataFile = dataFile;
  metadataFile += ".json";

  nlohmann::json metadata;
  metadata["num_events"] = numEvents;

  std::ofstream metadataStream(metadataFile);
  metadataStream << metadata.dump(2);
}

// Helper function to generate output filename
std::string generateOutputFilename(const std::string& prefix,
                                   std::size_t fileNum,
                                   HepMC3Util::Format format,
                                   HepMC3Util::Compression compression) {
  std::ostringstream filename;
  filename << prefix << "_" << std::setfill('0') << std::setw(6) << fileNum;

  if (format == HepMC3Util::Format::ascii) {
    filename << ".hepmc3";
  } else if (format == HepMC3Util::Format::root) {
    filename << ".root";
  }

  filename << HepMC3Util::compressionExtension(compression);

  return filename.str();
}

}  // namespace

HepMC3Util::NormalizeResult HepMC3Util::normalizeFiles(
    const std::vector<std::filesystem::path>& inputFiles,
    std::optional<std::filesystem::path> singleOutputPath,
    const std::filesystem::path& outputDir, const std::string& outputPrefix,
    std::optional<std::size_t> eventsPerFile,
    std::optional<std::size_t> maxEvents, Format format,
    Compression compression, int compressionLevel, bool verbose) {
  // Validate configuration
  if (inputFiles.empty()) {
    throw std::invalid_argument("No input files specified");
  }

  if (!singleOutputPath && !eventsPerFile.has_value()) {
    throw std::invalid_argument(
        "events-per-file must be > 0 in multi-file mode");
  }

  if (compressionLevel < 0 || compressionLevel > 19) {
    throw std::invalid_argument("compression-level must be 0-19");
  }

  if (format == Format::root && compression != Compression::none) {
    throw std::invalid_argument(
        "ROOT format does not support compression parameter");
  }

  NormalizeResult result;

  // Determine output directory
  std::filesystem::path actualOutputDir;
  if (singleOutputPath) {
    if (singleOutputPath->has_parent_path()) {
      actualOutputDir = singleOutputPath->parent_path();
    } else {
      actualOutputDir = ".";
    }
  } else {
    actualOutputDir = outputDir;
  }

  // Create output directory
  std::filesystem::create_directories(actualOutputDir);

  if (verbose) {
    std::cerr << "Configuration:\n";
    std::cerr << "  Input files: " << inputFiles.size() << "\n";
    if (singleOutputPath) {
      std::cerr << "  Single output file: " << *singleOutputPath << "\n";
      std::cerr << "  Format: " << format << "\n";
      std::cerr << "  Compression: " << compression << " (level "
                << compressionLevel << ")\n";
    } else {
      std::cerr << "  Output dir: " << actualOutputDir << "\n";
      std::cerr << "  Output prefix: " << outputPrefix << "\n";
      std::cerr << "  Events per file: "
                << (eventsPerFile.has_value()
                        ? std::to_string(eventsPerFile.value())
                        : "unset")
                << "\n";
      std::cerr << "  Format: " << format << "\n";
      std::cerr << "  Compression: " << compression << " (level "
                << compressionLevel << ")\n";
    }
    if (maxEvents.has_value()) {
      std::cerr << "  Max events to read: " << maxEvents.value() << "\n";
    }
    std::cerr << "\n";
  }

  // Processing state
  std::size_t globalEventIndex = 0;
  std::size_t eventsInCurrentFile = 0;
  std::size_t currentFileIndex = 0;
  std::unique_ptr<HepMC3::Writer> currentWriter;
  std::filesystem::path currentOutputPath;

  // Process each input file
  for (const auto& inputFile : inputFiles) {
    if (verbose) {
      std::cerr << "Reading " << inputFile << "...\n";
    }

    if (!std::filesystem::exists(inputFile)) {
      std::cerr << "WARNING: File not found: " << inputFile << "\n";
      continue;
    }

    // Track input file size
    result.totalInputSize += std::filesystem::file_size(inputFile);

    auto reader = HepMC3Util::deduceReader(inputFile.string());
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
      if (maxEvents > 0 && globalEventIndex >= maxEvents) {
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
      if (verbose &&
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

            std::size_t events =
                singleOutputPath ? globalEventIndex : eventsPerFile.value();
            writeMetadata(currentOutputPath, events);
            result.outputFiles.push_back(currentOutputPath);

            if (verbose) {
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
        if (singleOutputPath) {
          currentOutputPath = *singleOutputPath;
        } else {
          std::string filename = generateOutputFilename(
              outputPrefix, currentFileIndex, format, compression);
          currentOutputPath = actualOutputDir / filename;
        }
        currentWriter =
            HepMC3Util::createWriter(currentOutputPath, format, compression);
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
      if (!singleOutputPath && eventsInCurrentFile >= eventsPerFile) {
        currentWriter->close();

        // Get file size
        if (std::filesystem::exists(currentOutputPath)) {
          auto fileSize = std::filesystem::file_size(currentOutputPath);
          result.totalOutputSize += fileSize;

          writeMetadata(currentOutputPath, eventsInCurrentFile);
          result.outputFiles.push_back(currentOutputPath);

          if (verbose) {
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
        currentFileIndex++;
      }
    }

    reader->close();

    if (verbose) {
      // Clear progress line and print final count
      std::cerr << "  Read " << eventsReadFromFile << " events from "
                << inputFile << "                    \n";
    }

    // Check if we've reached the maximum number of events
    if (maxEvents.has_value() && globalEventIndex >= maxEvents.value()) {
      if (verbose) {
        std::cerr << "Reached maximum event limit (" << maxEvents.value()
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

      if (verbose) {
        double sizePerEvent = static_cast<double>(fileSize) /
                              static_cast<double>(eventsInCurrentFile);

        std::cerr << "  Wrote " << currentOutputPath << " ("
                  << eventsInCurrentFile << " events, " << formatSize(fileSize)
                  << ", "
                  << formatSize(static_cast<std::uintmax_t>(sizePerEvent))
                  << "/event)\n";
      }
    }
  }

  result.numEvents = globalEventIndex;

  // Print summary
  if (verbose) {
    std::cerr << "\nSummary:\n";
    if (singleOutputPath) {
      std::cerr << "  Processed " << globalEventIndex
                << " events into single file\n";
    } else {
      std::size_t totalFiles = (globalEventIndex + eventsPerFile.value() - 1) /
                               eventsPerFile.value();
      std::cerr << "  Processed " << globalEventIndex << " events into "
                << totalFiles << " file(s)\n";
    }

    if (globalEventIndex > 0) {
      double bytesPerEvent = static_cast<double>(result.totalOutputSize) /
                             static_cast<double>(globalEventIndex);
      std::cerr << "  Total input size:  " << formatSize(result.totalInputSize)
                << "\n";
      std::cerr << "  Total output size: " << formatSize(result.totalOutputSize)
                << " ("
                << formatSize(static_cast<std::uintmax_t>(bytesPerEvent))
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

namespace {
int eventGeneratorIndexImpl(const auto& obj) {
  return obj
      .template attribute<HepMC3::IntAttribute>(
          std::string{HepMC3Util::kEventGeneratorIndexAttribute})
      ->value();
}
}  // namespace

int HepMC3Util::eventGeneratorIndex(const HepMC3::GenParticle& particle) {
  return eventGeneratorIndexImpl(particle);
}

int HepMC3Util::eventGeneratorIndex(const HepMC3::GenVertex& vertex) {
  return eventGeneratorIndexImpl(vertex);
}

Acts::Vector4 HepMC3Util::convertPosition(const HepMC3::FourVector& vec) {
  using namespace Acts::UnitLiterals;
  return Acts::Vector4(vec.x() * 1_mm, vec.y() * 1_mm, vec.z() * 1_mm,
                       vec.t() * 1_mm);
}

std::vector<const HepMC3::GenVertex*> HepMC3Util::findHardScatterVertices(
    const HepMC3::GenEvent& event) {
  std::vector<const HepMC3::GenVertex*> vertices;

  constexpr double primaryVertexSpatialThreshold = 1 * Acts::UnitConstants::nm;

  for (const auto& vertex : event.vertices()) {
    // Convention is that idx=0 is hard-scatter
    if (HepMC3Util::eventGeneratorIndex(*vertex) == 0) {
      if (vertex->particles_in().empty() ||
          std::ranges::all_of(vertex->particles_in(), [](const auto& particle) {
            return particle->status() == HepMC3Util::kBeamParticleStatus;
          })) {
        auto it = std::ranges::find_if(vertices, [&](const auto& v) {
          return (HepMC3Util::convertPosition(v->position()) -
                  HepMC3Util::convertPosition(vertex->position()))
                     .template head<3>()
                     .cwiseAbs()
                     .maxCoeff() < primaryVertexSpatialThreshold;
        });

        if (it == vertices.end()) {
          vertices.push_back(vertex.get());
        }
      }
    }
  }

  return vertices;
}

}  // namespace ActsExamples
