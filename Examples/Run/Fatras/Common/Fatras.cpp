// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Fatras.hpp"

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>
#include <string>

#include <boost/program_options.hpp>

// helper functions are implemented in separate compilation units to reduce
// resource consumption during compilation

void addInputOptions(ActsExamples::Options::Description& desc);
void addSimulationOptions(ActsExamples::Options::Description& desc);
void addDigitizationOptions(ActsExamples::Options::Description& desc);
void setupInput(
    const ActsExamples::Options::Variables& variables,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers);
void setupSimulation(
    const ActsExamples::Options::Variables& variables,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry);
void setupDigitization(
    const ActsExamples::Options::Variables& variables,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry);

// fatras main function

int runFatras(int argc, char* argv[],
              std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  using namespace ActsExamples;

  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addOutputOptions(desc);
  addInputOptions(desc);
  // add general and detector-specific geometry options
  Options::addGeometryOptions(desc);
  detector->addOptions(desc);
  Options::addMaterialOptions(desc);
  // simulation also handles magnetic field
  addSimulationOptions(desc);
  addDigitizationOptions(desc);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic services
  auto randomNumbers =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));
  // make sure the output directory exists
  ActsExamples::ensureWritableDirectory(vars["output-dir"].as<std::string>());

  // setup sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));
  // setup detector geometry and material
  auto [trackingGeometry, contextDecorators] = Geometry::build(vars, *detector);
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }
  // setup other modules
  setupInput(vars, sequencer, randomNumbers);
  setupSimulation(vars, sequencer, randomNumbers, trackingGeometry);
  setupDigitization(vars, sequencer, randomNumbers, trackingGeometry);

  // run the simulation
  return sequencer.run();
}
