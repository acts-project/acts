// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "FatrasMain.hpp"

#include <boost/program_options.hpp>
#include <memory>

#include "ACTFW/Detector/IBaseDetector.hpp"
#include "ACTFW/Fatras/FatrasOptions.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Generators/ParticleSelector.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Options/ParticleGunOptions.hpp"
#include "ACTFW/Options/Pythia8Options.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "FatrasDigitizationBase.hpp"
#include "FatrasEvgenBase.hpp"
#include "FatrasSimulationBase.hpp"

int FW::fatrasMain(int argc, char* argv[],
                   std::shared_ptr<FW::IBaseDetector> detector) {
  using boost::program_options::value;

  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addParticleGunOptions(desc);
  FW::Options::addPythia8Options(desc);
  FW::Options::addRandomNumbersOptions(desc);
  FW::Options::addBFieldOptions(desc);
  FW::ParticleSelector::addOptions(desc);
  FW::Options::addFatrasOptions(desc);
  FW::Options::addOutputOptions(desc);
  desc.add_options()("evg-input-type",
                     value<std::string>()->default_value("pythia8"),
                     "Type of evgen input 'gun', 'pythia8'");
  // Add specific options for this geometry
  detector->addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  FW::Sequencer sequencer(FW::Options::readSequencerConfig(vm));

  // auto logLevel = FW::Options::readLogLevel(vm);

  // Create the random number engine
  auto randomNumberSvcCfg = FW::Options::readRandomNumbersConfig(vm);
  auto randomNumberSvc =
      std::make_shared<FW::RandomNumbers>(randomNumberSvcCfg);

  // The geometry, material and decoration
  auto geometry = FW::Geometry::build(vm, *detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // make sure the output directory exists
  FW::ensureWritableDirectory(vm["output-dir"].as<std::string>());

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
