// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Mille/MillePedeError.hpp"
#include "ActsPlugins/Mille/MillePedeSolver.hpp"

using namespace ActsPlugins::ActsToMille;

BOOST_AUTO_TEST_SUITE(MillePedeSolverBuiltin)

/// run the built-in silicon alignment test
BOOST_AUTO_TEST_CASE(BuiltInSiTracker) {
  MillePedeSolver::Config testSolverCfg;
  // this tells the solver we do not use a steering file
  testSolverCfg.steeringFile = "";
  testSolverCfg.workDir = "mpTest";
  testSolverCfg.redirectStdout = "mpStdout.txt";
  // this triggers the built-in test case
  testSolverCfg.extraOpts = {" -t=BRLF"};

  MillePedeSolver solver;
  auto res = solver.solve(testSolverCfg);
  // this one should be ok!
  BOOST_CHECK(res.ok());
  // the internal exit code should be zero
  BOOST_CHECK_EQUAL(res->exitCode, 0);
  // and the zero exit code should resolve to a "nominal exit"
  BOOST_CHECK_EQUAL(
      static_cast<int>(res->exitStatus),
      static_cast<int>(MillePedeSolver::mpExitStatus::nominalExit));
}

BOOST_AUTO_TEST_SUITE_END()
