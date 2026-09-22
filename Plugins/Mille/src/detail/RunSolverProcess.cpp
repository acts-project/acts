// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/detail/RunSolverProcess.hpp"

#include "ActsPlugins/Mille/detail/WrappedFileHandle.hpp"

#include <filesystem>

#include <boost/version.hpp>

// The mechanics of launching a child in Boost::Process vary
// between versions before and after 1.88. Make sure to
// pick up the one supported by the user's boost install.
#if BOOST_VERSION >= 108800
#include "SolverProcessCall_BoostV2.hpp"
#else
#include "SolverProcessCall_BoostV1.hpp"
#endif

ActsPlugins::ActsToMille::ChildProcessStatus
ActsPlugins::ActsToMille::runSolverProcess(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir,
    const std::filesystem::path& outputDest) {
  // work dir has to exist - caller's responsibility
  if (!std::filesystem::exists(runDir)) {
    return ActsPlugins::ActsToMille::ChildProcessStatus::FailedWorkDir;
  }

  WrappedFileHandle outputHandle;  // defaults to "do not redirect"

  if (!outputDest.empty()) {
    outputHandle = WrappedFileHandle(outputDest);
    if (!outputHandle.isRedirected()) {
      return ActsPlugins::ActsToMille::ChildProcessStatus::FailedRedirectStdout;
    }
  }
  return runChildProcessBoost(program, args, runDir, outputHandle);
}
