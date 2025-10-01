// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// @file HepMC3Normalize.cpp
/// @brief Standalone tool to normalize and chunk HepMC files

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <regex>
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

/// Detect compression type from file extension
std::string detectCompression(const std::filesystem::path& path) {
  std::string ext = path.extension().string();
  if (ext == ".gz") {
    return "zlib";
  } else if (ext == ".xz") {
    return "lzma";
  } else if (ext == ".bz2") {
    return "bzip2";
  } else if (ext == ".zst") {
    return "zstd";
  }
  return "none";
}

/// Detect format from file path
std::string detectFormat(const std::filesystem::path& path) {
  std::string pathStr = path.string();

  // Check for ROOT format: ends with .root
  if (std::regex_search(pathStr, std::regex(R"(\.root$)"))) {
    return "root";
  }

  // Check for HepMC3 ASCII format: ends with .hepmc3[.compression_ext] or
  // .hepmc[.compression_ext] Allowed compression extensions: .gz, .xz, .bz2,
  // .zst
  if (std::regex_search(pathStr,
                        std::regex(R"(\.hepmc3?(?:\.(?:gz|xz|bz2|zst))?$)"))) {
    return "ascii";
  }

  // If no pattern matches, throw an error
  throw std::invalid_argument(
      "Unable to detect format from filename: " + pathStr +
      "\nExpected patterns: *.hepmc3[.gz|.xz|.bz2|.zst], "
      "*.hepmc[.gz|.xz|.bz2|.zst], or *.root");
}

/// Generate output filename based on mode and parameters
std::string generateOutputFilename(const std::string& prefix,
                                   const std::string& format,
                                   const std::string& compression,
                                   std::size_t eventIndex) {
  // Multi-file mode: append event index
  std::string filename = std::format("{}_{:06d}", prefix, eventIndex);

  // Add format extension
  if (format == "ascii") {
    filename += ".hepmc3";
    filename += compressionExtension(compression);
  } else {
    filename += ".root";
  }

  return filename;
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

void printHelp(const po::options_description& desc, std::string_view command) {
  std::cerr << desc << "\n";
  std::cerr << "Examples:\n";
  std::cerr << "  # Multi-file mode with chunking\n";
  std::cerr << "  " << command
            << " -i input1.hepmc3.gz input2.hepmc3 -n 1000 -c zstd -l 9\n";
  std::cerr << "  " << command
            << " -i file.root -o normalized/ -p events -n 5000\n\n";
  std::cerr << "  # Single output mode (format/compression auto-detected)\n";
  std::cerr << "  " << command
            << " -i input1.hepmc3.gz input2.hepmc3 -S combined.hepmc3.zst\n";
  std::cerr << "  " << command
            << " -i input.hepmc3.gz -S output/events.hepmc3.gz\n";
}

int main(int argc, char** argv) {
  try {
    // Define command-line options
    po::options_description desc("HepMC3 File Normalizer\n\nOptions");
    // clang-format off
    desc.add_options()
      ("help,h", "Show this help message")
      ("list-compressions",
        "List available compression modes and exit")
      ("input,i",
        po::value<std::vector<std::string>>()->multitoken(),
        "Input HepMC files")
      ("single-output,S",
        po::value<std::string>(),
        "Write all events to a single output file. Format and compression are detected from filename.")
      ("output-dir,o",
        po::value<std::string>()->default_value("."),
        "Output directory (ignored with --single-output)")
      ("output-prefix,p",
        po::value<std::string>()->default_value("events"),
        "Output file prefix")
      ("events-per-file,n",
        po::value<std::size_t>()->default_value(10000),
        "Number of events per output file (ignored with --single-output)")
      ("max-events,m",
        po::value<std::size_t>()->default_value(0),
        "Maximum number of events to read (0 = all events)")
      ("compression,c",
        po::value<std::string>()->default_value("none"),
        "Compression type: none, zlib, lzma, bzip2, zstd (ignored with --single-output)")
      ("compression-level,l",
        po::value<int>()->default_value(6),
        "Compression level (higher = more compression)")
      ("format,f",
        po::value<std::string>()->default_value("ascii"),
       "Output format: ascii or root (ignored with --single-output)")
      ("json,j",
        "Write JSON output with list of created files to stdout")
      ("verbose,v", "Verbose output");
    // clang-format on

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help") > 0) {
      printHelp(desc, argv[0]);
      return 0;
    }

    if (vm.count("list-compressions") > 0) {
      std::cerr << "Available compression modes:\n";
      std::cerr << "  none       - No compression\n";
#ifdef HEPMC3_Z_SUPPORT
      std::cerr << "  zlib       - gzip compression (.gz)\n";
#endif
#ifdef HEPMC3_LZMA_SUPPORT
      std::cerr << "  lzma       - LZMA compression (.xz)\n";
#endif
#ifdef HEPMC3_BZ2_SUPPORT
      std::cerr << "  bzip2      - bzip2 compression (.bz2)\n";
#endif
#ifdef HEPMC3_ZSTD_SUPPORT
      std::cerr << "  zstd       - Zstandard compression (.zst)\n";
#endif
      std::cerr << "\nAvailable formats:\n";
      std::cerr << "  ascii      - HepMC3 ASCII format\n";
#ifdef ACTS_HEPMC3_ROOT_SUPPORT
      std::cerr << "  root       - ROOT format (has its own compression)\n";
#endif
      return 0;
    }

    po::notify(vm);

    if (vm.count("input") == 0) {
      std::cerr << "ERROR: No input files specified\n";
      printHelp(desc, argv[0]);
      return 1;
    }

    // Extract configuration
    auto inputFiles = vm["input"].as<std::vector<std::string>>();
    bool singleOutput = vm.count("single-output") > 0;
    std::filesystem::path singleOutputPath;
    std::filesystem::path outputDir;
    std::string outputPrefix;
    std::string compression;
    std::string format;

    if (singleOutput) {
      // Single output mode: detect format and compression from filename
      singleOutputPath = vm["single-output"].as<std::string>();
      compression = detectCompression(singleOutputPath);
      format = detectFormat(singleOutputPath);
      // Create parent directory if needed
      if (singleOutputPath.has_parent_path()) {
        outputDir = singleOutputPath.parent_path();
      } else {
        outputDir = ".";
      }
    } else {
      // Multi-file mode: use provided options
      outputDir = std::filesystem::path(vm["output-dir"].as<std::string>());
      outputPrefix = vm["output-prefix"].as<std::string>();
      compression = vm["compression"].as<std::string>();
      format = vm["format"].as<std::string>();
    }

    auto eventsPerFile = vm["events-per-file"].as<std::size_t>();
    auto maxEvents = vm["max-events"].as<std::size_t>();
    auto compressionLevel = vm["compression-level"].as<int>();
    bool verbose = vm.count("verbose") > 0;
    bool writeJson = vm.count("json") > 0;

    // Validate inputs
    if (inputFiles.empty()) {
      std::cerr << "ERROR: No input files specified\n";
      return 1;
    }

    if (!singleOutput && eventsPerFile == 0) {
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

    if (format == "root" && compression != "none") {
      std::cerr
          << "ERROR: ROOT format does not support compression parameter\n";
      return 1;
    }

    // Create output directory
    std::filesystem::create_directories(outputDir);

    if (verbose) {
      std::cerr << "Configuration:\n";
      std::cerr << "  Input files: " << inputFiles.size() << "\n";
      if (singleOutput) {
        std::cerr << "  Single output file: " << singleOutputPath << "\n";
        std::cerr << "  Format (detected): " << format << "\n";
        std::cerr << "  Compression (detected): " << compression << " (level "
                  << compressionLevel << ")\n";
      } else {
        std::cerr << "  Output dir: " << outputDir << "\n";
        std::cerr << "  Output prefix: " << outputPrefix << "\n";
        std::cerr << "  Events per file: " << eventsPerFile << "\n";
        std::cerr << "  Format: " << format << "\n";
        std::cerr << "  Compression: " << compression << " (level "
                  << compressionLevel << ")\n";
      }
      if (maxEvents > 0) {
        std::cerr << "  Max events to read: " << maxEvents << "\n";
      }
      std::cerr << "\n";
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

    // Timing tracking
    double totalWriteTime = 0.0;
    double totalReadTime = 0.0;

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
      constexpr std::size_t progressInterval = 100;

      while (!reader->failed()) {
        // Check if we've reached the maximum number of events
        if (maxEvents > 0 && globalEventIndex >= maxEvents) {
          break;
        }

        auto readStart = std::chrono::high_resolution_clock::now();
        reader->read_event(event);
        auto readEnd = std::chrono::high_resolution_clock::now();
        totalReadTime +=
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

            // Get file size and calculate compression ratio
            if (std::filesystem::exists(currentOutputPath)) {
              auto fileSize = std::filesystem::file_size(currentOutputPath);

              writeMetadata(currentOutputPath, eventsPerFile);
              outputFiles.push_back(currentOutputPath);

              if (verbose) {
                std::cerr << "  Wrote " << currentOutputPath << " ("
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

                double sizePerEvent = static_cast<double>(fileSize) /
                                      static_cast<double>(eventsPerFile);
                const char* perEventUnit = "B";
                if (sizePerEvent >= 1024.0) {
                  sizePerEvent /= 1024.0;
                  perEventUnit = "KiB";
                  if (sizePerEvent >= 1024.0) {
                    sizePerEvent /= 1024.0;
                    perEventUnit = "MiB";
                  }
                }

                std::cerr << std::fixed << std::setprecision(2) << size << " "
                          << unit << ", " << sizePerEvent << " " << perEventUnit
                          << "/event)\n";
              }
            }
          }

          // Generate output path based on mode
          if (singleOutput) {
            currentOutputPath = singleOutputPath;
          } else {
            std::string filename = generateOutputFilename(
                outputPrefix, format, compression, globalEventIndex);
            currentOutputPath = outputDir / filename;
          }
          currentWriter = createWriter(currentOutputPath, format, compression,
                                       compressionLevel);
        }

        // Set event number
        event.set_event_number(globalEventIndex);

        // Clear run info from events that are not the first in their output
        // file
        if (eventsInCurrentFile > 0) {
          event.set_run_info(nullptr);
        }

        auto writeStart = std::chrono::high_resolution_clock::now();
        currentWriter->write_event(event);
        auto writeEnd = std::chrono::high_resolution_clock::now();
        totalWriteTime +=
            std::chrono::duration<double>(writeEnd - writeStart).count();

        globalEventIndex++;
        eventsInCurrentFile++;
        eventsReadFromFile++;

        // Close file if chunk is complete (only in multi-file mode)
        if (!singleOutput && eventsInCurrentFile >= eventsPerFile) {
          currentWriter->close();

          // Get file size
          if (std::filesystem::exists(currentOutputPath)) {
            auto fileSize = std::filesystem::file_size(currentOutputPath);

            writeMetadata(currentOutputPath, eventsInCurrentFile);
            outputFiles.push_back(currentOutputPath);

            if (verbose) {
              std::cerr << "  Wrote " << currentOutputPath << " ("
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

              double sizePerEvent = static_cast<double>(fileSize) /
                                    static_cast<double>(eventsInCurrentFile);
              const char* perEventUnit = "B";
              if (sizePerEvent >= 1024.0) {
                sizePerEvent /= 1024.0;
                perEventUnit = "KiB";
                if (sizePerEvent >= 1024.0) {
                  sizePerEvent /= 1024.0;
                  perEventUnit = "MiB";
                }
              }

              std::cerr << std::format("{:.2f} {}, {:.2f} {}/event)\n", size,
                                       unit, sizePerEvent, perEventUnit);
            }
          }

          currentWriter.reset();
          eventsInCurrentFile = 0;
        }
      }

      reader->close();

      if (verbose) {
        // Clear progress line and print final count
        std::cerr << "  Read " << eventsReadFromFile << " events from "
                  << inputFile << "                    \n";
      }

      // Check if we've reached the maximum number of events
      if (maxEvents > 0 && globalEventIndex >= maxEvents) {
        if (verbose) {
          std::cerr << "Reached maximum event limit (" << maxEvents << ")\n";
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

        writeMetadata(currentOutputPath, eventsInCurrentFile);
        outputFiles.push_back(currentOutputPath);

        if (verbose) {
          std::cerr << "  Wrote " << currentOutputPath << " ("
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

          double sizePerEvent = static_cast<double>(fileSize) /
                                static_cast<double>(eventsInCurrentFile);
          const char* perEventUnit = "B";
          if (sizePerEvent >= 1024.0) {
            sizePerEvent /= 1024.0;
            perEventUnit = "KiB";
            if (sizePerEvent >= 1024.0) {
              sizePerEvent /= 1024.0;
              perEventUnit = "MiB";
            }
          }

          std::cerr << std::fixed << std::setprecision(2) << size << " " << unit
                    << ", " << sizePerEvent << " " << perEventUnit
                    << "/event)\n";
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
    std::cerr << "\nSummary:\n";
    if (singleOutput) {
      std::cerr << "  Processed " << globalEventIndex
                << " events into single file\n";
    } else {
      std::size_t totalFiles =
          (globalEventIndex + eventsPerFile - 1) / eventsPerFile;
      std::cerr << "  Processed " << globalEventIndex << " events into "
                << totalFiles << " file(s)\n";
    }

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

    if (globalEventIndex > 0) {
      double bytesPerEvent = static_cast<double>(totalOutputSize) /
                             static_cast<double>(globalEventIndex);
      std::cerr << "  Total input size:  " << formatSize(totalInputSize)
                << "\n";
      std::cerr << "  Total output size: " << formatSize(totalOutputSize)
                << " ("
                << formatSize(static_cast<std::uintmax_t>(bytesPerEvent))
                << "/event)\n";

      if (totalInputSize > 0) {
        double ratio = static_cast<double>(totalOutputSize) /
                       static_cast<double>(totalInputSize);
        std::cerr << "  Compression ratio: " << std::fixed
                  << std::setprecision(2) << (ratio * 100.0) << "%\n";
      }
    } else {
      std::cerr << "  Total input size:  " << formatSize(totalInputSize)
                << "\n";
      std::cerr << "  Total output size: " << formatSize(totalOutputSize)
                << "\n";
    }

    // Print timing information
    if (verbose && globalEventIndex > 0) {
      std::cerr << "\nTiming breakdown:\n";
      std::cerr << "  Reading events:  " << std::fixed << std::setprecision(3)
                << totalReadTime << " s ("
                << (totalReadTime / globalEventIndex * 1000.0)
                << " ms/event)\n";
      std::cerr << "  Writing events:  " << std::fixed << std::setprecision(3)
                << totalWriteTime << " s ("
                << (totalWriteTime / globalEventIndex * 1000.0)
                << " ms/event)\n";

      double totalProcessingTime = totalReadTime + totalWriteTime;
      std::cerr << "  Total processing: " << std::fixed << std::setprecision(3)
                << totalProcessingTime << " s\n";

      // Show percentage breakdown
      if (totalProcessingTime > 0) {
        std::cerr << "\nTime distribution:\n";
        std::cerr << "  Reading: " << std::fixed << std::setprecision(1)
                  << (totalReadTime / totalProcessingTime * 100.0) << "%\n";
        std::cerr << "  Writing: " << std::fixed << std::setprecision(1)
                  << (totalWriteTime / totalProcessingTime * 100.0) << "%\n";
      }
    }

    // Write JSON output if requested
    if (writeJson) {
      nlohmann::json j;
      j["num_events"] = globalEventIndex;
      j["num_files"] = outputFiles.size();
      j["files"] = nlohmann::json::array();

      for (const auto& file : outputFiles) {
        nlohmann::json fileInfo;
        // Use absolute path
        fileInfo["path"] = std::filesystem::absolute(file).string();
        if (std::filesystem::exists(file)) {
          fileInfo["size"] = std::filesystem::file_size(file);
        }
        j["files"].push_back(fileInfo);
      }

      std::cout << j.dump(2) << std::endl;
    }

    return 0;

  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
}
