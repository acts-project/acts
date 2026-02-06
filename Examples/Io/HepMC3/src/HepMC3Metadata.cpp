// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Metadata.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <fstream>

#include <nlohmann/json.hpp>

namespace ActsExamples::HepMC3Metadata {

static const std::string kEventCountKey = "num_events";

std::filesystem::path getSidecarPath(const std::filesystem::path& hepmc3File) {
  return hepmc3File.string() + ".json";
}

std::optional<HepMC3Metadata> readSidecar(
    const std::filesystem::path& hepmc3File) {
  auto sidecarPath = getSidecarPath(hepmc3File);

  if (!std::filesystem::exists(sidecarPath)) {
    return std::nullopt;
  }

  try {
    std::ifstream file(sidecarPath);
    if (!file.is_open()) {
      return std::nullopt;
    }

    nlohmann::json j;
    file >> j;

    if (!j.contains(kEventCountKey) ||
        !j[kEventCountKey].is_number_unsigned()) {
      return std::nullopt;
    }

    return HepMC3Metadata{.numEvents = j[kEventCountKey].get<std::size_t>()};
  } catch (...) {
    // Any error reading or parsing the file
    return std::nullopt;
  }
}

bool writeSidecar(const std::filesystem::path& hepmc3File,
                  const HepMC3Metadata& metadata, const Acts::Logger& logger) {
  auto sidecarPath = getSidecarPath(hepmc3File);
  ACTS_DEBUG("Writing HepMC3 sidecar metadata to " << sidecarPath);

  try {
    nlohmann::json j;
    j[kEventCountKey] = metadata.numEvents;

    std::ofstream file(sidecarPath);
    if (!file.is_open()) {
      ACTS_DEBUG("Failed to open sidecar file for writing: " << sidecarPath);
      return false;
    }

    file << j.dump(2);  // Pretty print with 2-space indent
    return file.good();
  } catch (...) {
    ACTS_DEBUG("Failed to write sidecar metadata file: " << sidecarPath);
    // Silent failure - directory might not be writable
    return false;
  }
}

}  // namespace ActsExamples::HepMC3Metadata
