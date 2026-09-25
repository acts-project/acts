// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/process.hpp>
#include <boost/process/environment.hpp>
#include <boost/process/process.hpp>
#include <boost/process/start_dir.hpp>
#include <boost/process/stdio.hpp>

namespace ActsPlugins::ActsToMille {
// With boost v1.88+, the v2 process API is default
ActsPlugins::ActsToMille::ChildProcessStatus runChildProcessBoost(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir,
    const WrappedFileHandle& outputHandle) {
  using namespace boost::process;

  // find the pede installation
  auto thePede = environment::find_executable(program);
  if (thePede.empty()) {
    return ActsPlugins::ActsToMille::ChildProcessStatus::ProgNotFound;
  }

  boost::asio::io_context io;

  process_stdio stdio{};
  if (outputHandle.isRedirected()) {
    stdio = process_stdio{{}, outputHandle(), outputHandle()};
  }
  // now run the fit
  process theProcess(io, thePede, args, process_start_dir(runDir.string()),
                     stdio);
  theProcess.wait();
  if (theProcess.exit_code() != 0) {
    return ChildProcessStatus::FailedRun;
  }
  return ChildProcessStatus::OK;
}
}  // namespace ActsPlugins::ActsToMille
