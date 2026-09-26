// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// Helper to define a "simple" interface for running pede
/// as a child process, without bringing in additional
/// external dependencies.

#include "Acts/Utilities/Logger.hpp"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace ActsPlugins {

enum class MpSolverStatus {
  OK = 0,
  ProgNotFound = 1,
  FailedRedirectStdout = 2,
  FailedWorkDir = 3,
  FailedRun = 4,
  UnknownError = 5
};

/// @brief attempt to run a program as a child process.
/// @param program: Program name (search on PATH)
/// @param args: command line args
/// @param runDir: Directory to run in. Caller is responsible for ensuring
/// validity.
/// @param logger: A logger instance
/// @param redirectOutput: if set, redirect pede output to file
/// of given name instead of printing to the host stdout
/// @return the call outcome as a status code
MpSolverStatus runSolverProcess(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir, const Acts::Logger& logger,
    const std::optional<std::filesystem::path>& redirectOutput = std::nullopt);

}  // namespace ActsPlugins
