// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Utilities/Units.hpp>
#include <cstdlib>
#include <memory>

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Generators/FlattenEvent.hpp"
#include "ACTFW/Io/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Io/Root/RootParticleWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Options/ParticleGunOptions.hpp"
#include "ACTFW/Printers/PrintParticles.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using namespace FW;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addParticleGunOptions(desc);
  Options::addOutputOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  Acts::Logging::Level logLevel = Options::readLogLevel(vm);

  // basic services
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // event generation w/ particle gun
  EventGenerator::Config evgenCfg = Options::readParticleGunOptions(vm);
  evgenCfg.output = "event";
  evgenCfg.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgenCfg, logLevel));

  // flatten event to just particles
  FlattenEvent::Config flatten;
  flatten.inputEvent = evgenCfg.output;
  flatten.outputParticles = "particles";
  sequencer.addAlgorithm(std::make_shared<FlattenEvent>(flatten, logLevel));

  // print generated particles
  PrintParticles::Config printCfg;
  printCfg.inputParticles = flatten.outputParticles;
  sequencer.addAlgorithm(std::make_shared<PrintParticles>(printCfg, logLevel));

  // different output modes
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  if (vm["output-csv"].as<bool>()) {
    CsvParticleWriter::Config csvWriterCfg;
    csvWriterCfg.inputParticles = flatten.outputParticles;
    csvWriterCfg.outputDir = outputDir;
    csvWriterCfg.outputStem = "particles";
    sequencer.addWriter(
        std::make_shared<CsvParticleWriter>(csvWriterCfg, logLevel));
  }
  if (vm["output-root"].as<bool>()) {
    RootParticleWriter::Config rootWriterCfg;
    rootWriterCfg.inputParticles = flatten.outputParticles;
    rootWriterCfg.filePath = joinPaths(outputDir, "particles.root");
    sequencer.addWriter(
        std::make_shared<RootParticleWriter>(rootWriterCfg, logLevel));
  }

  return sequencer.run();
}
