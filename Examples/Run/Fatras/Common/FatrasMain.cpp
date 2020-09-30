// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "FatrasMain.hpp"

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Fatras/FatrasOptions.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/Options/Pythia8Options.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

#include <boost/program_options.hpp>

#include "FatrasDigitizationBase.hpp"
#include "FatrasEvgenBase.hpp"
#include "FatrasSimulationBase.hpp"

int ActsExamples::fatrasMain(
    int argc, char* argv[],
    std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  using boost::program_options::value;

  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addParticleGunOptions(desc);
  ActsExamples::Options::addPythia8Options(desc);
  ActsExamples::Options::addRandomNumbersOptions(desc);
  ActsExamples::Options::addBFieldOptions(desc);
  ActsExamples::ParticleSelector::addOptions(desc);
  ActsExamples::Options::addFatrasOptions(desc);
  ActsExamples::Options::addDigitizationOptions(desc);
  ActsExamples::Options::addOutputOptions(desc);
  desc.add_options()("evg-input-type",
                     value<std::string>()->default_value("pythia8"),
                     "Type of evgen input 'gun', 'pythia8'");
  // Add specific options for this geometry
  detector->addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // Create the random number engine
  auto randomNumberSvcCfg = ActsExamples::Options::readRandomNumbersConfig(vm);
  auto randomNumberSvc =
      std::make_shared<ActsExamples::RandomNumbers>(randomNumberSvcCfg);

  // The geometry, material and decoration
  auto geometry = ActsExamples::Geometry::build(vm, *detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // make sure the output directory exists
  ActsExamples::ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // (A) EVGEN
  // Setup the evgen input to the simulation
  setupEvgenInput(vm, sequencer, randomNumberSvc);

  // (B) SIMULATION
  // Setup the simulation
  setupSimulation(vm, sequencer, randomNumberSvc, tGeometry);

  // (C) DIGITIZATION
  // Setup the digitization
  setupDigitization(vm, sequencer, randomNumberSvc, tGeometry);

  return sequencer.run();
}
