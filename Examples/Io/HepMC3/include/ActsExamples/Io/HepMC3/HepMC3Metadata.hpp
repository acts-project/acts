// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>

namespace Acts {
class Logger;
}

namespace ActsExamples::HepMC3Metadata {

/// Get the path to the sidecar metadata file for a given HepMC3 file.
/// Example: "file.hepmc3" -> "file.hepmc3.meta"
///          "file.root" -> "file.root.meta"
std::filesystem::path getSidecarPath(const std::filesystem::path& hepmc3File);

/// Read the event count from a sidecar metadata file.
/// Returns std::nullopt if the file doesn't exist, is invalid, or can't be
/// read.
std::optional<std::size_t> readSidecar(const std::filesystem::path& hepmc3File);

/// Write the event count to a sidecar metadata file.
/// Returns true if successful, false if the write failed (e.g., directory not
/// writable). Failures are silent and do not throw exceptions.
bool writeSidecar(const std::filesystem::path& hepmc3File,
                  std::size_t eventCount, const Acts::Logger& logger);

}  // namespace ActsExamples::HepMC3Metadata
