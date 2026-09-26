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

enum class ChildProcessStatus {
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
/// @param outputDest: if not empty, redirect output to a file instead
/// of printing to the host stdout
/// @return the call outcome as a status code
ChildProcessStatus runSolverProcess(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir,
    const std::filesystem::path& outputDest = "");

}  // namespace ActsPlugins::ActsToMille
