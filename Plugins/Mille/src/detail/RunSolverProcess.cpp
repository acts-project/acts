// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "RunSolverProcess.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <filesystem>

#include "WrappedFileHandle.hpp"

// The mechanics of launching a child in Boost::Process vary
// between versions before and after 1.88. Make sure to
// pick up the one supported by the user's boost install.
#ifdef ACTS_MILLE_USE_BOOST_PROCESS_V2
#include "SolverProcessCall_BoostV2.hpp"
#else
#include "SolverProcessCall_BoostV1.hpp"
#endif

namespace ActsPlugins::detail {
MpSolverStatus runSolverProcess(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir, const Acts::Logger& logger,
    const std::optional<std::filesystem::path>& redirectOutput) {
  // work dir has to exist - caller's responsibility
  if (!std::filesystem::exists(runDir)) {
    ACTS_ERROR("Run directory '" << runDir
                                 << "' does not exist, will not run solver");
    return MpSolverStatus::FailedWorkDir;
  }

  WrappedFileHandle outputHandle;  // defaults to "do not redirect"

  if (redirectOutput.has_value()) {
    outputHandle = WrappedFileHandle(*redirectOutput);
    if (!outputHandle.isRedirected()) {
      ACTS_ERROR("Failed to redirect output to '" << *redirectOutput << "'");
      return MpSolverStatus::FailedRedirectStdout;
    }
  }
  return runChildProcessBoost(program, args, runDir, outputHandle, logger);
}
}  // namespace ActsPlugins::detail
