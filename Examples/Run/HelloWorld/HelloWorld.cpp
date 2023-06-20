// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file
/// @brief An example tools that shows the sequencer functionality

#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"

#include <cstdlib>
#include <memory>

#include "HelloLoggerAlgorithm.hpp"
#include "HelloRandomAlgorithm.hpp"
#include "HelloWhiteBoardAlgorithm.hpp"

int main(int argc, char* argv[]) {
  // setup options
  // every component should have an associated option setup function
  // that should be called here.
  auto opt = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(opt);
  ActsExamples::Options::addRandomNumbersOptions(opt);
  // parse options from command line flags
  auto vm = ActsExamples::Options::parse(opt, argc, argv);
  // an empty variables map indicates an error
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // extract some common options
  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // setup basic tools shared among algorithms
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
      ActsExamples::Options::readRandomNumbersConfig(vm));

  // setup the sequencer first w/ config derived from options
  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // add HelloWorld algorithm that does nothing
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::HelloLoggerAlgorithm>(logLevel));

  // add HelloRandom algorithm that uses RandomNumbers to generate some
  // random numbers from various distributions.
  ActsExamples::HelloRandomAlgorithm::Config rndCfg;
  rndCfg.randomNumbers = rnd;
  rndCfg.gaussParameters = {{0., 2.5}};
  rndCfg.uniformParameters = {{-1.23, 4.25}};
  rndCfg.gammaParameters = {{1., 1.}};
  rndCfg.drawsPerEvent = 5000;
  rndCfg.output = "random_data";
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::HelloRandomAlgorithm>(rndCfg, logLevel));

  // add HelloWhiteBoardAlgorithm the reads/writes data from/to the event store
  ActsExamples::HelloWhiteBoardAlgorithm::Config wbCfg;
  // use data from previous algorithm as input
  wbCfg.input = rndCfg.output;
  wbCfg.output = "copied_data";
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::HelloWhiteBoardAlgorithm>(wbCfg,
                                                               logLevel));

  // Run all configured algorithms and return the appropriate status.
  return sequencer.run();
}
