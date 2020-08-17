// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "FatrasEvgenBase.hpp"

#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Generators/FlattenEvent.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/Options/Pythia8Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

void ActsExamples::setupEvgenInput(
    const ActsExamples::Options::Variables& vm,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumberSvc) {
  // Read the standard options
  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // Add requested event generator
  auto evgenInput = vm["evg-input-type"].template as<std::string>();
  if (evgenInput == "gun") {
    auto evgCfg = ActsExamples::Options::readParticleGunOptions(vm);
    evgCfg.output = "event_generated";
    evgCfg.randomNumbers = randomNumberSvc;
    sequencer.addReader(
        std::make_shared<ActsExamples::EventGenerator>(evgCfg, logLevel));

  } else if (evgenInput == "pythia8") {
    auto evgCfg = ActsExamples::Options::readPythia8Options(vm, logLevel);
    evgCfg.output = "event_generated";
    evgCfg.randomNumbers = randomNumberSvc;
    sequencer.addReader(
        std::make_shared<ActsExamples::EventGenerator>(evgCfg, logLevel));

  } else {
    throw std::runtime_error("unknown event generator input: " + evgenInput);
  }

  // Convert to particles for the writers and subsequent algorithms
  ActsExamples::FlattenEvent::Config flatten;
  flatten.inputEvent = "event_generated";
  flatten.outputParticles = "particles_generated";
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::FlattenEvent>(flatten, logLevel));

  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();

  // Write generated particles as CSV files
  if (vm["output-csv"].template as<bool>()) {
    ActsExamples::CsvParticleWriter::Config pWriterCsvConfig;
    pWriterCsvConfig.inputParticles = flatten.outputParticles;
    pWriterCsvConfig.outputDir = outputDir;
    pWriterCsvConfig.outputStem = "particles_generated";
    sequencer.addWriter(std::make_shared<ActsExamples::CsvParticleWriter>(
        pWriterCsvConfig, logLevel));
  }

  // Write generated particles as ROOT file
  if (vm["output-root"].template as<bool>()) {
    // Write particles as ROOT TTree
    ActsExamples::RootParticleWriter::Config pWriterRootConfig;
    pWriterRootConfig.inputParticles = flatten.outputParticles;
    pWriterRootConfig.filePath =
        ActsExamples::joinPaths(outputDir, "particles.root");
    sequencer.addWriter(std::make_shared<ActsExamples::RootParticleWriter>(
        pWriterRootConfig, logLevel));
  }
}
