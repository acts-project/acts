// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Mille/detail/runChildProcess.hpp"

#include <fstream>

using namespace ActsPlugins::ActsToMille;

BOOST_AUTO_TEST_SUITE(ChildProcessTests)

/// catch a missing steering file
BOOST_AUTO_TEST_CASE(MissingProg) {
  childProcessStatus status =
      runChildProcess("tHiSsur3lyD0esn.otExist.exe", {});
  BOOST_CHECK(status == childProcessStatus::progNotFound);
}

/// catch a missing steering file
BOOST_AUTO_TEST_CASE(WrongWD) {
  childProcessStatus status =
      runChildProcess("echo", {"Hello World"}, "/non/existing/directory");
  BOOST_CHECK(status == childProcessStatus::failedWorkDir);
}

/// successful call
BOOST_AUTO_TEST_CASE(GoodCall) {
  childProcessStatus status = runChildProcess("echo", {"Hello World"});
  BOOST_CHECK(status == childProcessStatus::ok);
}

/// redirect stdout
BOOST_AUTO_TEST_CASE(Redirect) {
  const std::string testMessage = "Hello ACTS!";
  childProcessStatus status =
      runChildProcess("echo", {testMessage}, "", "teststdout.txt");
  BOOST_CHECK(status == childProcessStatus::ok);
  std::ifstream in("teststdout.txt");
  BOOST_CHECK(in.is_open());
  std::string read = "";
  BOOST_CHECK(std::getline(in, read).good());
  BOOST_CHECK(read == testMessage);
}

BOOST_AUTO_TEST_SUITE_END()
