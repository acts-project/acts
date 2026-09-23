// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Diagnostics.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
#include "ActsPlugins/Json/TrackingGeometryMaterialJsonConverter.hpp"
#include "ActsPlugins/Json/detail/JsonIo.hpp"

#include <charconv>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
constexpr std::string_view usage =
    "Usage: ActsMaterialMapMigrate INPUT OUTPUT [OPTIONS]\n"
    "Migrate deprecated surface material maps to the version 1 format.\n"
    "Input: JSON or CBOR, optionally zstd compressed (detected from content).\n"
    "Output: .json, .cbor, .json.zst or .cbor.zst.\n"
    "Volume material is not supported by version 1.\n\n"
    "Options:\n"
    "  --material-fraction-bits N  Retained fraction bits, 0-23 (default 23)\n"
    "  --compression-level N       zstd level (default 9)\n"
    "  --indentation N             JSON indentation (default 4)\n"
    "  --help                      Show this help\n";

template <typename T>
T number(std::string_view text) {
  T result{};
  const auto [end, error] =
      std::from_chars(text.data(), text.data() + text.size(), result);
  if (error != std::errc{} || end != text.data() + text.size()) {
    throw std::invalid_argument("Invalid numeric option: " + std::string(text));
  }
  return result;
}
}  // namespace

int main(int argc, char* argv[]) {
  for (int i = 1; i < argc; ++i) {
    if (std::string_view(argv[i]) == "--help") {
      std::cout << usage;
      return 0;
    }
  }
  if (argc < 3) {
    std::cerr << usage;
    return 1;
  }
  try {
    Acts::TrackingGeometryMaterialJsonConverter::Options options;
    for (int i = 3; i < argc; ++i) {
      const std::string_view option = argv[i];
      if (++i == argc) {
        throw std::invalid_argument("Missing value for " + std::string(option));
      }
      if (option == "--material-fraction-bits") {
        options.materialFractionBits = number<unsigned int>(argv[i]);
      } else if (option == "--compression-level") {
        options.compressionLevel = number<int>(argv[i]);
      } else if (option == "--indentation") {
        options.indentation = number<unsigned int>(argv[i]);
      } else {
        throw std::invalid_argument("Unknown option: " + std::string(option));
      }
    }
    const std::filesystem::path input = argv[1];
    const std::filesystem::path output = argv[2];
    if (std::filesystem::exists(output) &&
        std::filesystem::equivalent(input, output)) {
      throw std::invalid_argument("Input and output must be different files");
    }
    const auto document = Acts::detail::readJsonFile(input);
    // Check the envelope before calling the legacy reader, which uses [].
    if (!document.contains("Surfaces") || !document.contains("Volumes")) {
      throw std::invalid_argument(
          "Expected a deprecated material map with Surfaces and Volumes");
    }
    // Migration intentionally reads the deprecated format.
    ACTS_PUSH_IGNORE_DEPRECATED()
    Acts::MaterialMapJsonConverter legacy({}, Acts::Logging::WARNING);
    const auto material = legacy.jsonToMaterialMaps(document);
    ACTS_POP_IGNORE_DEPRECATED()
    Acts::TrackingGeometryMaterialJsonConverter().toFile(material, output,
                                                         options);
    std::cout << "Migrated "
              << material.surfaceMaterials.size() +
                     material.keyedSurfaces.size()
              << " surface assignments to " << output << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Material migration failed: " << error.what() << '\n';
    return 1;
  }
}
