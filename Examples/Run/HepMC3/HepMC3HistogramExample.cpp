// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/HepMC/HepMCProcessExtractor.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3Reader.hpp"
#include "ActsExamples/Io/NuclearInteractions/RootNuclearInteractionParametersWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/HepMC3Options.hpp"

///
/// Straight forward example of reading a HepMC3 file.
///
int main(int argc, char** argv) {
  // Declare the supported program options.
  // Setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addInputOptions(desc);
  ActsExamples::Options::addHepMC3ReaderOptions(desc);

  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // Create the reader
  auto hepMC3ReaderConfig = ActsExamples::Options::readHepMC3ReaderOptions(vm);
  hepMC3ReaderConfig.outputEvents = "hepmc-events";

  ActsExamples::HepMCProcessExtractor::Config extractionConfig;
  extractionConfig.inputEvents = hepMC3ReaderConfig.outputEvents;
  extractionConfig.extractionProcess = "Inelastic";

  ActsExamples::RootNuclearInteractionParametersWriter::Config writerCfg;
  writerCfg.inputSimulationProcesses =
      extractionConfig.outputSimulationProcesses;

  // Add to the sequencer
  sequencer.addReader(std::make_shared<ActsExamples::HepMC3AsciiReader>(
      hepMC3ReaderConfig, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::HepMCProcessExtractor>(
      std::move(extractionConfig), logLevel));
  sequencer.addWriter(
      std::make_shared<ActsExamples::RootNuclearInteractionParametersWriter>(
          writerCfg, logLevel));

  // Run
  return sequencer.run();
}
