// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// @file HepMC3Normalize.cpp
/// @brief Standalone CLI tool to normalize and chunk HepMC files

#include "ActsExamples/Io/HepMC3/HepMC3Normalizer.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

namespace po = boost::program_options;

using namespace ActsExamples;

/// Detect compression type from file extension
HepMC3Util::Compression detectCompression(const std::filesystem::path& path) {
  std::string ext = path.extension().string();
  if (ext == ".gz") {
    return HepMC3Util::Compression::zlib;
  } else if (ext == ".xz") {
    return HepMC3Util::Compression::lzma;
  } else if (ext == ".bz2") {
    return HepMC3Util::Compression::bzip2;
  } else if (ext == ".zst") {
    return HepMC3Util::Compression::zstd;
  }
  return HepMC3Util::Compression::none;
}

/// Convert string to Compression enum
HepMC3Util::Compression stringToCompression(const std::string& str) {
  if (str == "none") {
    return HepMC3Util::Compression::none;
  } else if (str == "zlib") {
    return HepMC3Util::Compression::zlib;
  } else if (str == "lzma") {
    return HepMC3Util::Compression::lzma;
  } else if (str == "bzip2") {
    return HepMC3Util::Compression::bzip2;
  } else if (str == "zstd") {
    return HepMC3Util::Compression::zstd;
  }
  throw std::invalid_argument("Unknown compression type: " + str);
}

/// Convert string to Format enum
HepMC3Util::Format stringToFormat(const std::string& str) {
  if (str == "ascii") {
    return HepMC3Util::Format::ascii;
  } else if (str == "root") {
    return HepMC3Util::Format::root;
  }
  throw std::invalid_argument("Unknown format: " + str);
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
      for (auto comp : HepMC3Util::availableCompressionModes()) {
        std::cerr << "  " << comp << "\n";
      }
      std::cerr << "\nAvailable formats:\n";
      for (auto fmt : HepMC3Util::availableFormats()) {
        std::cerr << "  " << fmt << "\n";
      }
      return 0;
    }

    po::notify(vm);

    if (vm.count("input") == 0) {
      std::cerr << "ERROR: No input files specified\n";
      printHelp(desc, argv[0]);
      return 1;
    }

    // Build configuration
    HepMC3Normalizer::Config cfg;

    // Input files
    auto inputStrings = vm["input"].as<std::vector<std::string>>();
    for (const auto& str : inputStrings) {
      cfg.inputFiles.emplace_back(str);
    }

    // Single output or multi-file mode
    bool singleOutput = vm.count("single-output") > 0;
    if (singleOutput) {
      std::filesystem::path singleOutputPath =
          vm["single-output"].as<std::string>();
      cfg.singleOutputPath = singleOutputPath;
      cfg.compression = detectCompression(singleOutputPath);
      cfg.format = HepMC3Util::formatFromFilename(singleOutputPath.string());
    } else {
      cfg.outputDir = vm["output-dir"].as<std::string>();
      cfg.outputPrefix = vm["output-prefix"].as<std::string>();
      cfg.eventsPerFile = vm["events-per-file"].as<std::size_t>();
      cfg.compression = stringToCompression(vm["compression"].as<std::string>());
      cfg.format = stringToFormat(vm["format"].as<std::string>());
    }

    cfg.maxEvents = vm["max-events"].as<std::size_t>();
    cfg.compressionLevel = vm["compression-level"].as<int>();
    cfg.verbose = vm.count("verbose") > 0;

    bool writeJson = vm.count("json") > 0;

    // Run normalization
    HepMC3Normalizer normalizer(cfg);
    auto result = normalizer.normalize();

    // Write JSON output if requested
    if (writeJson) {
      nlohmann::json j;
      j["num_events"] = result.numEvents;
      j["num_files"] = result.outputFiles.size();
      j["files"] = nlohmann::json::array();

      for (const auto& file : result.outputFiles) {
        nlohmann::json fileInfo;
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
