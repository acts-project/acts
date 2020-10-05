// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file An example utility to read and print particles CSV files.
///
/// This examples shows how to use the framework CSV I/O and printer algorithms.

#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <memory>

#include <boost/program_options.hpp>

using namespace ActsExamples;
using boost::program_options::value;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions("Read and print particles CSVs");
  auto opts = desc.add_options();
  Options::addSequencerOptions(desc);
  opts("input-dir", value<std::string>()->default_value(""), "");
  opts("input-stem", value<std::string>()->default_value("particles"), "");
  ParticleSelector::addOptions(desc);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vars));

  // read options
  auto logLevel = Options::readLogLevel(vars);
  auto inputDir = vars["input-dir"].as<std::string>();
  auto inputStem = vars["input-stem"].as<std::string>();

  // read particles
  CsvParticleReader::Config readParticlesCfg;
  readParticlesCfg.inputDir = inputDir;
  readParticlesCfg.inputStem = inputStem;
  readParticlesCfg.outputParticles = "particles";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(readParticlesCfg, logLevel));

  // pre-select particles
  auto selectParticlesCfg = ParticleSelector::readConfig(vars);
  selectParticlesCfg.inputParticles = readParticlesCfg.outputParticles;
  selectParticlesCfg.outputParticles = "particles_selected";
  sequencer.addAlgorithm(
      std::make_shared<ParticleSelector>(selectParticlesCfg, logLevel));

  // print selected particles
  ParticlesPrinter::Config printParticlesCfg;
  printParticlesCfg.inputParticles = selectParticlesCfg.outputParticles;
  sequencer.addAlgorithm(
      std::make_shared<ParticlesPrinter>(printParticlesCfg, logLevel));

  return sequencer.run();
}
