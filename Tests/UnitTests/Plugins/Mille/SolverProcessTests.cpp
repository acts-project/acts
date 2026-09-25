// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Mille/detail/RunSolverProcess.hpp"

#include <filesystem>
#include <fstream>

using namespace ActsPlugins;

BOOST_AUTO_TEST_SUITE(ChildProcessTests)
auto logger = Acts::getDefaultLogger("MillePedeSolver", Acts::Logging::INFO);
/// catch a missing steering file
BOOST_AUTO_TEST_CASE(MissingProg) {
// temporarily disable log failure threshold
// as we intentionally test an error condition
#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  auto level = Acts::Logging::getFailureThreshold();
  Acts::Logging::setFailureThreshold(Acts::Logging::MAX);
#endif
  MpSolverStatus status =
      runSolverProcess("tHiSsur3lyD0esn.otExist.exe", {},
                       std::filesystem::current_path(), *logger);
  BOOST_CHECK(status == MpSolverStatus::ProgNotFound);
#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  Acts::Logging::setFailureThreshold(level);
#endif
}

/// catch an invalid work dir
BOOST_AUTO_TEST_CASE(WrongWD) {
// temporarily disable log failure threshold
// as we intentionally test an error condition
#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  auto level = Acts::Logging::getFailureThreshold();
  Acts::Logging::setFailureThreshold(Acts::Logging::MAX);
#endif
  MpSolverStatus status = runSolverProcess(
      "echo", {"Hello World"}, "o/hNo/Invalid/WorkDirectory", *logger);
  BOOST_CHECK(status == MpSolverStatus::FailedWorkDir);
#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
  Acts::Logging::setFailureThreshold(level);
#endif
}

/// successful call
BOOST_AUTO_TEST_CASE(GoodCall) {
  MpSolverStatus status = runSolverProcess(
      "echo", {"Hello World"}, std::filesystem::current_path(), *logger);
  BOOST_CHECK(status == MpSolverStatus::OK);
}

/// redirect stdout
BOOST_AUTO_TEST_CASE(Redirect) {
  const std::string testMessage = "Hello ACTS!";
  MpSolverStatus status =
      runSolverProcess("echo", {testMessage}, std::filesystem::current_path(),
                       *logger, "teststdout.txt");
  BOOST_CHECK(status == MpSolverStatus::OK);
  std::ifstream in("teststdout.txt");
  BOOST_CHECK(in.is_open());
  std::string read = "";
  BOOST_CHECK(std::getline(in, read).good());
  BOOST_CHECK(read == testMessage);
}

BOOST_AUTO_TEST_SUITE_END()
