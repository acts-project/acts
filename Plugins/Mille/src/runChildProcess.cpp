// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/detail/runChildProcess.hpp"

#include <array>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace ActsPlugins::ActsToMille;

childProcessStatus ActsPlugins::ActsToMille::runChildProcess(
    const std::string& program, const std::vector<std::string>& args,
    const std::filesystem::path& runDir,
    const std::filesystem::path& output_dest) {
  // support for redirecting pede stdout
  int stdout_fd = -1;
  if (!output_dest.empty()) {
    stdout_fd = ::open(output_dest.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (stdout_fd == -1) {
      std::cerr << " Failed to redirect output to `" << output_dest << "`"
                << std::endl;
      return childProcessStatus::failedRedirectStdout;
    }
  }

  /// setup a pipe to pass error information
  /// from the child to the parent
  std::array<int, 2> errCodePipe{0, 0};

  if (::pipe(errCodePipe.data()) == -1) {
    std::cerr << " Failed to open pipe for error codes" << std::endl;
    return childProcessStatus::failedFork;
  }
  // Make the write end close automatically if the call in our child succeeds
  if (::fcntl(errCodePipe[1], F_SETFD, FD_CLOEXEC) == -1) {
    ::close(errCodePipe[0]);
    ::close(errCodePipe[1]);
    std::cerr << " Failed to open pipe for error codes" << std::endl;
    return childProcessStatus::failedFork;
  }
  // fork child process for pede
  const pid_t pid = ::fork();

  // fork failed
  if (pid == -1) {
    if (stdout_fd != -1) {
      ::close(stdout_fd);
    }
    ::close(errCodePipe[0]);
    ::close(errCodePipe[1]);
    std::cerr << " Failed to fork child process with error " << errno
              << std::endl;
    return childProcessStatus::failedFork;
  }

  /// the pede child process
  if (pid == 0) {
    ::close(errCodePipe[0]);
    // redirect I/O if configured by the user
    if (stdout_fd != -1) {
      if (::dup2(stdout_fd, STDOUT_FILENO) == -1 ||
          ::dup2(stdout_fd, STDERR_FILENO) == -1) {
        std::cerr << " Failed to redirect I/O " << std::endl;
      }
      ::close(stdout_fd);
      stdout_fd = -1;
      // intentionally continue, with output to standard streams
    }

    // change run directory
    if (!runDir.empty()) {
      std::cout << " try to CD into " << runDir << std::endl;
      if (::chdir(runDir.c_str()) == -1) {
        childProcessStatus fail = childProcessStatus::failedWorkDir;
        (void)::write(errCodePipe[1], &fail, sizeof(fail));
        ::_exit(999);
      }
    }

    // build CLI arg list
    std::vector<char*> argv;
    argv.reserve(args.size() + 2);

    argv.push_back(const_cast<char*>(program.data()));

    for (const auto& arg : args) {
      argv.push_back(const_cast<char*>(arg.c_str()));
    }

    argv.push_back(nullptr);

    if (::execvp(argv[0], argv.data()) == -1) {
      if (errno == ENOENT || errno == ENOEXEC) {
        childProcessStatus fail = childProcessStatus::progNotFound;
        (void)::write(errCodePipe[1], &fail, sizeof(fail));
      } else {
        childProcessStatus fail = childProcessStatus::failedRun;
        (void)::write(errCodePipe[1], &fail, sizeof(fail));
      }
      ::_exit(999);
    }
  }

  // Parent

  ::close(errCodePipe[1]);
  childProcessStatus childStat = childProcessStatus::ok;
  const ssize_t nRead = ::read(errCodePipe[0], &childStat, sizeof(childStat));
  ::close(errCodePipe[0]);

  // wait for pede to finish
  int status = 0;

  if (::waitpid(pid, &status, 0) == -1) {
    std::cerr << " Failed to run child " << std::endl;
    return childProcessStatus::failedRun;
  }
  // close the targetr file for redirected output, if it exists
  if (stdout_fd != -1) {
    ::close(stdout_fd);
  }
  if (nRead == sizeof(childStat)) {
    return childStat;
  }

  if (WIFEXITED(status)) {
    return childProcessStatus::ok;
  }

  if (WIFSIGNALED(status)) {
    std::cerr << "process terminated by signal " +
                     std::to_string(WTERMSIG(status))
              << std::endl;
    return childProcessStatus::caughtSignal;
  }

  std::cerr << "Failed to run child" << std::endl;
  return childProcessStatus::unknownError;
}
