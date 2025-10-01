// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// @file HepMC3Normalize.cpp
/// @brief Standalone tool to normalize and chunk HepMC files

#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderFactory.h>
#include <HepMC3/Writer.h>
#include <HepMC3/WriterAscii.h>
#include <boost/program_options.hpp>

#ifdef HEPMC3_USE_COMPRESSION
#include <HepMC3/CompressedIO.h>
#endif

#ifdef ACTS_HEPMC3_ROOT_SUPPORT
#include <HepMC3/WriterRootTree.h>
#endif

#include <nlohmann/json.hpp>

namespace po = boost::program_options;

/// Helper to get file extension for compression type
std::string compressionExtension(const std::string& compression) {
  if (compression == "none") {
    return "";
  } else if (compression == "zlib") {
    return ".gz";
  } else if (compression == "lzma") {
    return ".xz";
  } else if (compression == "bzip2") {
    return ".bz2";
  } else if (compression == "zstd") {
    return ".zst";
  }
  return "";
}

#ifdef HEPMC3_USE_COMPRESSION
/// Convert string to bxz::Compression enum
bxz::Compression stringToCompression(const std::string& compression) {
  if (compression == "none") {
    return bxz::Compression::plaintext;
  }
#ifdef HEPMC3_Z_SUPPORT
  else if (compression == "zlib") {
    return bxz::Compression::z;
  }
#endif
#ifdef HEPMC3_LZMA_SUPPORT
  else if (compression == "lzma") {
    return bxz::Compression::lzma;
  }
#endif
#ifdef HEPMC3_BZ2_SUPPORT
  else if (compression == "bzip2") {
    return bxz::Compression::bz2;
  }
#endif
#ifdef HEPMC3_ZSTD_SUPPORT
  else if (compression == "zstd") {
    return bxz::Compression::zstd;
  }
#endif
  throw std::invalid_argument("Unsupported compression type: " + compression);
}
#endif

/// Write metadata sidecar file
void writeMetadata(const std::filesystem::path& hepmcFile,
                   std::size_t eventCount) {
  auto metaPath = hepmcFile.string() + ".meta";

  try {
    nlohmann::json j;
    j["eventCount"] = eventCount;

    std::ofstream out(metaPath);
    if (out.is_open()) {
      out << j.dump(2);
    }
  } catch (...) {
    // Silent failure - not critical
  }
}

/// Copy event with full attribute preservation
void copyEvent(const HepMC3::GenEvent& src, HepMC3::GenEvent& dst) {
  // Copy basic properties
  dst.set_event_number(src.event_number());
  dst.set_units(src.momentum_unit(), src.length_unit());

  // Copy event-level attributes
  for (const auto& attr : src.attribute_names()) {
    auto value = src.attribute_as_string(attr);
    dst.add_attribute(attr, std::make_shared<HepMC3::StringAttribute>(value));
  }

  // Copy particles (keep indexed for topology reconstruction)
  std::vector<std::shared_ptr<HepMC3::GenParticle>> particles;
  particles.reserve(src.particles().size());

  for (const auto& srcParticle : src.particles()) {
    auto particle = std::make_shared<HepMC3::GenParticle>();
    particle->set_momentum(srcParticle->momentum());
    particle->set_generated_mass(srcParticle->generated_mass());
    particle->set_pid(srcParticle->pid());
    particle->set_status(srcParticle->status());

    // Copy particle attributes
    for (const auto& attr : srcParticle->attribute_names()) {
      auto value = srcParticle->attribute_as_string(attr);
      particle->add_attribute(attr,
                              std::make_shared<HepMC3::StringAttribute>(value));
    }

    particles.push_back(particle);
    dst.add_particle(particle);
  }

  // Copy vertices with topology
  for (const auto& srcVertex : src.vertices()) {
    auto vertex = std::make_shared<HepMC3::GenVertex>(srcVertex->position());
    vertex->set_status(srcVertex->status());

    // Copy vertex attributes
    for (const auto& attr : srcVertex->attribute_names()) {
      auto value = srcVertex->attribute_as_string(attr);
      vertex->add_attribute(attr,
                            std::make_shared<HepMC3::StringAttribute>(value));
    }

    dst.add_vertex(vertex);

    // Reconnect particle topology using indexed particle array
    for (const auto& srcParticle : srcVertex->particles_in()) {
      vertex->add_particle_in(particles.at(srcParticle->id() - 1));
    }
    for (const auto& srcParticle : srcVertex->particles_out()) {
      vertex->add_particle_out(particles.at(srcParticle->id() - 1));
    }
  }
}

/// Create HepMC3 writer with optional compression
std::unique_ptr<HepMC3::Writer> createWriter(const std::filesystem::path& path,
                                             const std::string& format,
                                             const std::string& compression,
                                             int compressionLevel) {
  // ROOT format
  if (format == "root") {
    if (compression != "none") {
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
  bxz::Compression comp = stringToCompression(compression);

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
  if (compression != "none") {
    throw std::runtime_error("Compression support not enabled in this build");
  }
  return std::make_unique<HepMC3::WriterAscii>(path.string());
#endif
}

int main(int argc, char** argv) {
  try {
    // Define command-line options
    po::options_description desc("HepMC3 File Normalizer\n\nOptions");
    desc.add_options()("help,h", "Show this help message")(
        "list-compressions", "List available compression modes and exit")(
        "input,i", po::value<std::vector<std::string>>(),
        "Input HepMC files (can be specified multiple times)")(
        "output-dir,o", po::value<std::string>()->default_value("."),
        "Output directory")("output-prefix,p",
                            po::value<std::string>()->default_value("events"),
                            "Output file prefix")(
        "events-per-file,n", po::value<std::size_t>()->default_value(10000),
        "Number of events per output file")(
        "compression,c", po::value<std::string>()->default_value("none"),
        "Compression type: none, zlib, lzma, bzip2, zstd")(
        "compression-level,l", po::value<int>()->default_value(6),
        "Compression level (0-9, higher = more compression)")(
        "format,f", po::value<std::string>()->default_value("ascii"),
        "Output format: ascii or root")("verbose,v", "Verbose output");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      std::cout << "Examples:\n";
      std::cout << "  " << argv[0]
                << " -i input1.hepmc3.gz -i input2.hepmc3 -n 1000 -c zstd -l "
                   "9\n";
      std::cout << "  " << argv[0]
                << " -i file.root -o normalized/ -p events -n 5000\n";
      return 0;
    }

    if (vm.count("list-compressions")) {
      std::cout << "Available compression modes:\n";
      std::cout << "  none       - No compression\n";
#ifdef HEPMC3_Z_SUPPORT
      std::cout << "  zlib       - gzip compression (.gz)\n";
#endif
#ifdef HEPMC3_LZMA_SUPPORT
      std::cout << "  lzma       - LZMA compression (.xz)\n";
#endif
#ifdef HEPMC3_BZ2_SUPPORT
      std::cout << "  bzip2      - bzip2 compression (.bz2)\n";
#endif
#ifdef HEPMC3_ZSTD_SUPPORT
      std::cout << "  zstd       - Zstandard compression (.zst)\n";
#endif
      std::cout << "\nAvailable formats:\n";
      std::cout << "  ascii      - HepMC3 ASCII format\n";
#ifdef ACTS_HEPMC3_ROOT_SUPPORT
      std::cout << "  root       - ROOT format (has its own compression)\n";
#endif
      return 0;
    }

    po::notify(vm);

    // Extract configuration
    auto inputFiles = vm["input"].as<std::vector<std::string>>();
    auto outputDir = std::filesystem::path(vm["output-dir"].as<std::string>());
    auto outputPrefix = vm["output-prefix"].as<std::string>();
    auto eventsPerFile = vm["events-per-file"].as<std::size_t>();
    auto compression = vm["compression"].as<std::string>();
    auto compressionLevel = vm["compression-level"].as<int>();
    auto format = vm["format"].as<std::string>();
    bool verbose = vm.count("verbose") > 0;

    // Validate inputs
    if (inputFiles.empty()) {
      std::cerr << "ERROR: No input files specified\n";
      return 1;
    }

    if (eventsPerFile == 0) {
      std::cerr << "ERROR: events-per-file must be > 0\n";
      return 1;
    }

    if (compressionLevel < 0 || compressionLevel > 19) {
      std::cerr << "ERROR: compression-level must be 0-9\n";
      return 1;
    }

    if (format != "ascii" && format != "root") {
      std::cerr << "ERROR: format must be 'ascii' or 'root'\n";
      return 1;
    }

    // Create output directory
    std::filesystem::create_directories(outputDir);

    if (verbose) {
      std::cout << "Configuration:\n";
      std::cout << "  Input files: " << inputFiles.size() << "\n";
      std::cout << "  Output dir: " << outputDir << "\n";
      std::cout << "  Events per file: " << eventsPerFile << "\n";
      std::cout << "  Format: " << format << "\n";
      std::cout << "  Compression: " << compression << " (level "
                << compressionLevel << ")\n\n";
    }

    // Processing state
    std::size_t globalEventIndex = 0;
    std::size_t eventsInCurrentFile = 0;
    std::unique_ptr<HepMC3::Writer> currentWriter;
    std::filesystem::path currentOutputPath;

    // File size tracking
    std::uintmax_t totalInputSize = 0;
    std::uintmax_t totalOutputSize = 0;
    std::vector<std::filesystem::path> outputFiles;

    // Process each input file
    for (const auto& inputFile : inputFiles) {
      if (verbose) {
        std::cout << "Reading " << inputFile << "...\n";
      }

      if (!std::filesystem::exists(inputFile)) {
        std::cerr << "WARNING: File not found: " << inputFile << "\n";
        continue;
      }

      // Track input file size
      totalInputSize += std::filesystem::file_size(inputFile);

      auto reader = HepMC3::deduce_reader(inputFile);
      if (!reader) {
        std::cerr << "ERROR: Failed to open " << inputFile << "\n";
        continue;
      }

      HepMC3::GenEvent event;
      std::size_t eventsReadFromFile = 0;
      std::size_t lastProgressUpdate = 0;
      // Update every 1000 events
      constexpr std::size_t progressInterval = 1000;

      while (!reader->failed()) {
        reader->read_event(event);
        if (reader->failed()) {
          break;
        }

        // Show progress
        if (verbose &&
            (eventsReadFromFile - lastProgressUpdate) >= progressInterval) {
          std::cout << "    Progress: " << eventsReadFromFile
                    << " events read...\r" << std::flush;
          lastProgressUpdate = eventsReadFromFile;
        }

        // Create new output file if needed
        if (eventsInCurrentFile == 0) {
          // Close previous file
          if (currentWriter) {
            currentWriter->close();

            // Get file size and calculate compression ratio
            if (std::filesystem::exists(currentOutputPath)) {
              auto fileSize = std::filesystem::file_size(currentOutputPath);

              writeMetadata(currentOutputPath, eventsPerFile);
              outputFiles.push_back(currentOutputPath);

              if (verbose) {
                std::cout << "  Wrote " << currentOutputPath << " ("
                          << eventsPerFile << " events, ";

                // Format file size
                double size = static_cast<double>(fileSize);
                const char* unit = "B";
                if (size >= 1024.0) {
                  size /= 1024.0;
                  unit = "KiB";
                  if (size >= 1024.0) {
                    size /= 1024.0;
                    unit = "MiB";
                    if (size >= 1024.0) {
                      size /= 1024.0;
                      unit = "GiB";
                    }
                  }
                }

                double sizePerEvent = static_cast<double>(fileSize) / static_cast<double>(eventsPerFile);
                const char* perEventUnit = "B";
                if (sizePerEvent >= 1024.0) {
                  sizePerEvent /= 1024.0;
                  perEventUnit = "KiB";
                  if (sizePerEvent >= 1024.0) {
                    sizePerEvent /= 1024.0;
                    perEventUnit = "MiB";
                  }
                }

                std::cout << std::fixed << std::setprecision(2) << size << " "
                          << unit << ", " << sizePerEvent << " " << perEventUnit << "/event)\n";
              }
            }
          }

          // Generate output filename: events_000000.hepmc3.zst
          std::ostringstream filename;
          filename << outputPrefix << "_" << std::setw(6) << std::setfill('0')
                   << globalEventIndex;

          if (format == "ascii") {
            filename << ".hepmc3";
            filename << compressionExtension(compression);
          } else {
            filename << ".root";
          }

          currentOutputPath = outputDir / filename.str();
          currentWriter = createWriter(currentOutputPath, format, compression,
                                       compressionLevel);
        }

        // Copy and write event
        HepMC3::GenEvent outEvent;
        copyEvent(event, outEvent);
        outEvent.set_event_number(globalEventIndex);
        currentWriter->write_event(outEvent);

        globalEventIndex++;
        eventsInCurrentFile++;
        eventsReadFromFile++;

        // Close file if chunk is complete
        if (eventsInCurrentFile >= eventsPerFile) {
          currentWriter->close();

          // Get file size
          if (std::filesystem::exists(currentOutputPath)) {
            auto fileSize = std::filesystem::file_size(currentOutputPath);

            writeMetadata(currentOutputPath, eventsInCurrentFile);
            outputFiles.push_back(currentOutputPath);

            if (verbose) {
              std::cout << "  Wrote " << currentOutputPath << " ("
                        << eventsInCurrentFile << " events, ";

              // Format file size
              double size = static_cast<double>(fileSize);
              const char* unit = "B";
              if (size >= 1024.0) {
                size /= 1024.0;
                unit = "KiB";
                if (size >= 1024.0) {
                  size /= 1024.0;
                  unit = "MiB";
                  if (size >= 1024.0) {
                    size /= 1024.0;
                    unit = "GiB";
                  }
                }
              }

              double sizePerEvent = static_cast<double>(fileSize) / static_cast<double>(eventsInCurrentFile);
              const char* perEventUnit = "B";
              if (sizePerEvent >= 1024.0) {
                sizePerEvent /= 1024.0;
                perEventUnit = "KiB";
                if (sizePerEvent >= 1024.0) {
                  sizePerEvent /= 1024.0;
                  perEventUnit = "MiB";
                }
              }

              std::cout << std::format("{:.2f} {}, {:.2f} {}/event)\n",
                                       size, unit, sizePerEvent, perEventUnit);
            }
          }

          currentWriter.reset();
          eventsInCurrentFile = 0;
        }
      }

      reader->close();

      if (verbose) {
        // Clear progress line and print final count
        std::cout << "  Read " << eventsReadFromFile << " events from "
                  << inputFile << "                    \n";
      }
    }

    // Close final file
    if (currentWriter && eventsInCurrentFile > 0) {
      currentWriter->close();

      // Get file size
      if (std::filesystem::exists(currentOutputPath)) {
        auto fileSize = std::filesystem::file_size(currentOutputPath);

        writeMetadata(currentOutputPath, eventsInCurrentFile);
        outputFiles.push_back(currentOutputPath);

        if (verbose) {
          std::cout << "  Wrote " << currentOutputPath << " ("
                    << eventsInCurrentFile << " events, ";

          // Format file size
          double size = static_cast<double>(fileSize);
          const char* unit = "B";
          if (size >= 1024.0) {
            size /= 1024.0;
            unit = "KiB";
            if (size >= 1024.0) {
              size /= 1024.0;
              unit = "MiB";
              if (size >= 1024.0) {
                size /= 1024.0;
                unit = "GiB";
              }
            }
          }

          double sizePerEvent = static_cast<double>(fileSize) / static_cast<double>(eventsInCurrentFile);
          const char* perEventUnit = "B";
          if (sizePerEvent >= 1024.0) {
            sizePerEvent /= 1024.0;
            perEventUnit = "KiB";
            if (sizePerEvent >= 1024.0) {
              sizePerEvent /= 1024.0;
              perEventUnit = "MiB";
            }
          }

          std::cout << std::fixed << std::setprecision(2) << size << " " << unit
                    << ", " << sizePerEvent << " " << perEventUnit << "/event)\n";
        }
      }
    }

    // Calculate total output size
    for (const auto& outFile : outputFiles) {
      if (std::filesystem::exists(outFile)) {
        totalOutputSize += std::filesystem::file_size(outFile);
      }
    }

    // Print summary
    std::size_t totalFiles =
        (globalEventIndex + eventsPerFile - 1) / eventsPerFile;
    std::cout << "\nSummary:\n";
    std::cout << "  Processed " << globalEventIndex << " events into "
              << totalFiles << " file(s)\n";

    // Format file sizes
    auto formatSize = [](std::uintmax_t bytes) -> std::string {
      const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
      int unitIndex = 0;
      double size = static_cast<double>(bytes);

      while (size >= 1024.0 && unitIndex < 4) {
        size /= 1024.0;
        unitIndex++;
      }

      std::ostringstream oss;
      oss << std::fixed << std::setprecision(2) << size << " "
          << units[unitIndex];
      return oss.str();
    };

    std::cout << "  Total input size:  " << formatSize(totalInputSize) << " ("
              << totalInputSize << " bytes)\n";
    std::cout << "  Total output size: " << formatSize(totalOutputSize) << " ("
              << totalOutputSize << " bytes)\n";

    if (totalInputSize > 0) {
      double ratio = static_cast<double>(totalOutputSize) /
                     static_cast<double>(totalInputSize);
      std::cout << "  Compression ratio: " << std::fixed << std::setprecision(2)
                << (ratio * 100.0) << "%\n";
    }

    return 0;

  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
