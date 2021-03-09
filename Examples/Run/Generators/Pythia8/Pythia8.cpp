// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/Pythia8Options.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Definitions/Units.hpp>

#include <cstdlib>
#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addPythia8Options(desc);
  Options::addOutputOptions(desc, OutputFormat::Csv | OutputFormat::Root);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = Options::readLogLevel(vm);
  Sequencer::Config sequencerCfg = Options::readSequencerConfig(vm);
  Sequencer sequencer(sequencerCfg);

  // basic services
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // event generation w/ internal pythia8 instance
  EventGenerator::Config evgen = Options::readPythia8Options(vm, logLevel);
  evgen.outputParticles = "particles";
  evgen.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgen, logLevel));

  // print generated particles
  if ((logLevel == Acts::Logging::VERBOSE) or
      (logLevel == Acts::Logging::DEBUG)) {
    ParticlesPrinter::Config print;
    print.inputParticles = evgen.outputParticles;
    sequencer.addAlgorithm(std::make_shared<ParticlesPrinter>(print, logLevel));
  }

  // different output modes
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  if (vm["output-csv"].as<bool>()) {
    CsvParticleWriter::Config csvWriter;
    csvWriter.inputParticles = evgen.outputParticles;
    csvWriter.outputDir = outputDir;
    csvWriter.outputStem = "particles";
    sequencer.addWriter(
        std::make_shared<CsvParticleWriter>(csvWriter, logLevel));
  }
  if (vm["output-root"].as<bool>()) {
    RootParticleWriter::Config rootWriter;
    rootWriter.inputParticles = evgen.outputParticles;
    rootWriter.filePath = joinPaths(outputDir, "particles.root");
    sequencer.addWriter(
        std::make_shared<RootParticleWriter>(rootWriter, logLevel));
  }

  return sequencer.run();
}
