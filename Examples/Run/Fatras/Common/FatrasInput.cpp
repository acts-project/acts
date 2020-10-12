// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <boost/program_options.hpp>

#include "FatrasInternal.hpp"

void addInputOptions(ActsExamples::Options::Description& desc) {
  using namespace ActsExamples;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("input-dir", value<std::string>(),
      "Read particle input from CSV files in the given directory. If not "
      "given, particles are generated with the particle gun");
  ActsExamples::Options::addParticleGunOptions(desc);
  ActsExamples::ParticleSelector::addOptions(desc);
}

void setupInput(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers) {
  using namespace ActsExamples;

  // Read the standard options
  auto logLevel = ActsExamples::Options::readLogLevel(vars);

  if (not vars["input-dir"].empty()) {
    // read particle input from csv file
    CsvParticleReader::Config readParticles;
    readParticles.inputDir = vars["input-dir"].as<std::string>();
    readParticles.inputStem = "particles";
    readParticles.outputParticles = "particles_input";
    sequencer.addReader(
        std::make_shared<CsvParticleReader>(readParticles, logLevel));

  } else {
    // generate particle input from a particle gun
    auto gen = Options::readParticleGunOptions(vars);
    gen.outputParticles = "particles_input";
    gen.randomNumbers = randomNumbers;
    sequencer.addReader(std::make_shared<EventGenerator>(gen, logLevel));
  }

  // add additional particle selection
  auto select = ActsExamples::ParticleSelector::readConfig(vars);
  select.inputParticles = "particles_input";
  select.outputParticles = kFatrasCollectionParticles;
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::ParticleSelector>(select, logLevel));
}
