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

#include <filesystem>
#include <string>
#include <vector>

namespace ActsPlugins::ActsToMille {

enum class childProcessStatus {
  ok = 0,
  progNotFound = 1,
  failedRedirectStdout = 2,
  failedWorkDir = 3,
  failedFork = 4,
  failedRun = 5,
  caughtSignal = 6,
  unknownError = 7
};

/// @brief attempt to run a program as a child process.
/// @param program: Program name (search on PATH)
/// @param args: command line args
/// @param path: if not empty, redirect output to a file instead
/// of printing to the host stdout
/// @param runDir: If not empty, will run in the specified directory.
/// Should already exist (caller is responsible).
/// @return a pair with the call outcome and the return code (if any, else -1).
childProcessStatus runChildProcess(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir = "",
    const std::filesystem::path& output_dest = "");

}  // namespace ActsPlugins::ActsToMille
