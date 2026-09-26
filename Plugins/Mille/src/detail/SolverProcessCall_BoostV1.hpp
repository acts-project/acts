// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/process.hpp>

namespace ActsPlugins::ActsToMille {
ActsPlugins::ActsToMille::ChildProcessStatus runChildProcessBoost(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir,
    const WrappedFileHandle& outputHandle) {
  namespace bp = boost::process;

  // find the pede installation
  auto thePede = bp::search_path(program);
  if (thePede.empty()) {
    return ActsPlugins::ActsToMille::ChildProcessStatus::ProgNotFound;
  }

  bp::child theProcess;
  if (outputHandle.isRedirected()) {
    theProcess = bp::child(thePede, args, bp::start_dir = runDir.string(),
                           (bp::std_out & bp::std_err) > outputHandle());
  } else {
    theProcess = bp::child(thePede, args, bp::start_dir = runDir.string());
  }

  theProcess.wait();
  if (theProcess.exit_code() != 0) {
    return ChildProcessStatus::FailedRun;
  }
  return ChildProcessStatus::OK;
}

}  // namespace ActsPlugins::ActsToMille
