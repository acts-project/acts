// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/detail/runChildProcess.hpp"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <boost/version.hpp>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#if BOOST_VERSION >= 108800

#include <boost/process/environment.hpp>
#include <boost/process/process.hpp>
#include <boost/process/start_dir.hpp>
#include <boost/process/stdio.hpp>

#else
#include <boost/process.hpp>
#endif

using namespace ActsPlugins::ActsToMille;

namespace {

/// @brief helper to wrap a file handle
class wrappedFileHandle {
 public:
  explicit wrappedFileHandle(const std::filesystem::path& outf = "") {
    if (!outf.empty()) {
      m_handle = open(outf.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    }
  }
  ~wrappedFileHandle() {
    if (m_handle >= 0) {
      close(m_handle);
    }
  }
  wrappedFileHandle(const wrappedFileHandle&) = delete;
  wrappedFileHandle& operator=(const wrappedFileHandle&) = delete;

  wrappedFileHandle(wrappedFileHandle&& other) noexcept
      : m_handle(std::exchange(other.m_handle, -1)) {}

  wrappedFileHandle& operator=(wrappedFileHandle&& other) noexcept {
    if (this != &other) {
      if (m_handle >= 0) {
        ::close(m_handle);
      }
      m_handle = std::exchange(other.m_handle, -1);
    }

    return *this;
  }
  int operator()() const { return m_handle; }
  bool isRedirected() const { return m_handle != -1; }

 private:
  int m_handle = -1;
};

#if BOOST_VERSION >= 108800
// With boost v1.88+, the v2 process API is default
ActsPlugins::ActsToMille::childProcessStatus runChildProcessBoost(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir,
    const wrappedFileHandle& outputHandle) {
  using namespace boost::process;

  // find the pede installation
  auto thePede = environment::find_executable(program);
  if (thePede.empty()) {
    return ActsPlugins::ActsToMille::childProcessStatus::progNotFound;
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
    std::cout << theProcess.exit_code() << std::endl;
    return childProcessStatus::failedRun;
  }
  return childProcessStatus::ok;
}

// earlier boost versions use the v1 API.
#else

#include <boost/process.hpp>

ActsPlugins::ActsToMille::childProcessStatus runChildProcessBoost(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir,
    const wrappedFileHandle& outputHandle) {
  using namespace boost::process;

  // find the pede installation
  auto thePede = search_path(program);
  if (thePede.empty()) {
    return ActsPlugins::ActsToMille::childProcessStatus::progNotFound;
  }

  std::unique_ptr<bp::child> theProcess = nullptr;
  if (outputHandle.isRedirected) {
    theProcess = std::make_unique<bp::child>(
        thePede, args, bp::start_dir = workDir.string(),
        bp::std_out > outputHandle(), bp::std_err > outputHandle());
  } else {
      theProcess = std::make_unique<bp::child>( thePede,
        args,
        bp::start_dir = workDir.string()
  };

  theProcess->wait();
  if (theProcess->exit_code() != 0) {
    return childProcessStatus::failedRun;
  }
  return childProcessStatus::ok;
}

#endif

}  // namespace

ActsPlugins::ActsToMille::childProcessStatus
ActsPlugins::ActsToMille::runChildProcess(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir,
    const std::filesystem::path& output_dest) {
  // determine where the user wishes to run
  std::filesystem::path workDir = std::filesystem::current_path();

  if (!runDir.empty()) {
    workDir = runDir;
  }
  if (!std::filesystem::exists(workDir)) {
    std::error_code e;
    std::filesystem::create_directories(workDir, e);
    if (e) {
      return ActsPlugins::ActsToMille::childProcessStatus::failedWorkDir;
    }
  };

  wrappedFileHandle outputHandle;  // defaults to "do not redirect"

  if (!output_dest.empty()) {
    outputHandle = wrappedFileHandle(output_dest);
    if (outputHandle() < 0) {
      return ActsPlugins::ActsToMille::childProcessStatus::failedRedirectStdout;
    }
  }

  return runChildProcessBoost(program, args, workDir, outputHandle);
}
