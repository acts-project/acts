// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Mille/detail/RunSolverProcess.hpp"

#include <filesystem>
#include <fstream>

using namespace ActsPlugins::ActsToMille;

BOOST_AUTO_TEST_SUITE(ChildProcessTests)

/// catch a missing steering file
BOOST_AUTO_TEST_CASE(MissingProg) {
  ChildProcessStatus status = runSolverProcess(
      "tHiSsur3lyD0esn.otExist.exe", {}, std::filesystem::current_path());
  BOOST_CHECK(status == ChildProcessStatus::ProgNotFound);
}

/// catch an invalid work dir
BOOST_AUTO_TEST_CASE(WrongWD) {
  ChildProcessStatus status =
      runSolverProcess("echo", {"Hello World"}, "o/hNo/Invalid/WorkDirectory");
  BOOST_CHECK(status == ChildProcessStatus::FailedWorkDir);
}

/// successful call
BOOST_AUTO_TEST_CASE(GoodCall) {
  ChildProcessStatus status = runSolverProcess("echo", {"Hello World"},
                                               std::filesystem::current_path());
  BOOST_CHECK(status == ChildProcessStatus::OK);
}

/// redirect stdout
BOOST_AUTO_TEST_CASE(Redirect) {
  const std::string testMessage = "Hello ACTS!";
  ChildProcessStatus status = runSolverProcess(
      "echo", {testMessage}, std::filesystem::current_path(), "teststdout.txt");
  BOOST_CHECK(status == ChildProcessStatus::OK);
  std::ifstream in("teststdout.txt");
  BOOST_CHECK(in.is_open());
  std::string read = "";
  BOOST_CHECK(std::getline(in, read).good());
  BOOST_CHECK(read == testMessage);
}

BOOST_AUTO_TEST_SUITE_END()
