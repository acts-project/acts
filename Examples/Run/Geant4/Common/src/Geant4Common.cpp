// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Common.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"
#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"
#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"
#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSteppingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"
#include "ActsExamples/Geant4/SimParticleTranslation.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/Simulation/CommonSimulation.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>
#include <string>

#include <FTFP_BERT.hh>
#include <G4RunManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <boost/program_options.hpp>

namespace ActsExamples {

void setupMaterialRecording(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<DetectorConstructionFactory> detectorConstructionFactory) {
  auto g4loglevel =
      Acts::Logging::Level(vars["g4-loglevel"].as<unsigned int>());
  size_t seed = vars["g4-seed"].as<size_t>();

  Geant4MaterialRecording::Config g4Cfg;

  g4Cfg.detectorConstructionFactory = std::move(detectorConstructionFactory);
  g4Cfg.randomNumbers = std::make_shared<ActsExamples::RandomNumbers>(
      ActsExamples::RandomNumbers::Config{seed});
  g4Cfg.inputParticles = Simulation::kParticlesInitial;
  g4Cfg.outputMaterialTracks = Simulation::kMaterialTracks;
  g4Cfg.excludeMaterials = {"Air", "Vacuum"};

  sequencer.addAlgorithm(
      std::make_shared<Geant4MaterialRecording>(g4Cfg, g4loglevel));
}

void setupGeant4Simulation(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<DetectorConstructionFactory> detectorConstructionFactory,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField) {
  auto g4loglevel =
      Acts::Logging::Level(vars["g4-loglevel"].as<unsigned int>());
  size_t seed = vars["g4-seed"].as<size_t>();

  Geant4Simulation::Config g4Cfg;

  g4Cfg.detectorConstructionFactory = std::move(detectorConstructionFactory);
  g4Cfg.randomNumbers = std::make_shared<ActsExamples::RandomNumbers>(
      ActsExamples::RandomNumbers::Config{seed});
  g4Cfg.inputParticles = Simulation::kParticlesSelection;
  g4Cfg.outputSimHits = Simulation::kSimHits;
  g4Cfg.outputParticlesInitial = Simulation::kParticlesInitial;
  g4Cfg.outputParticlesFinal = Simulation::kParticlesFinal;
  g4Cfg.trackingGeometry = std::move(trackingGeometry);
  g4Cfg.magneticField = std::move(magneticField);

  sequencer.addAlgorithm(std::make_shared<Geant4Simulation>(g4Cfg, g4loglevel));
}

int runMaterialRecording(
    const ActsExamples::Options::Variables& vars,
    std::shared_ptr<DetectorConstructionFactory> detectorConstructionFactory) {
  // Set up the sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));

  // Output level and log level
  auto logLevel = Options::readLogLevel(vars);

  std::string outputDir =
      ensureWritableDirectory(vars["output-dir"].as<std::string>());

  // Basic services
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));

  // Event generation w/ particle gun
  EventGenerator::Config evgen = Options::readParticleGunOptions(vars);
  evgen.outputParticles = Simulation::kParticlesInitial;
  evgen.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgen, logLevel));

  // Set up the Geant4 Simulation
  setupMaterialRecording(vars, sequencer,
                         std::move(detectorConstructionFactory));

  // setup the output writing
  if (vars["output-root"].as<bool>()) {
    // Write the propagation steps as ROOT TTree
    RootMaterialTrackWriter::Config materialTrackWriter;
    materialTrackWriter.prePostStep = true;
    materialTrackWriter.recalculateTotals = true;
    materialTrackWriter.collection = Simulation::kMaterialTracks;
    materialTrackWriter.filePath = joinPaths(
        outputDir,
        "geant4_" + std::string(Simulation::kMaterialTracks) + ".root");
    sequencer.addWriter(std::make_shared<RootMaterialTrackWriter>(
        materialTrackWriter, logLevel));
  }

  auto result = sequencer.run();
  return result;
}

int runGeant4Simulation(
    const ActsExamples::Options::Variables& vars,
    std::shared_ptr<DetectorConstructionFactory> detectorConstructionFactory,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
  // Basic services
  auto randomNumbers =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));

  // Set up the sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));

  // Set up the magnetic field
  auto magneticField = ActsExamples::Options::readMagneticField(vars);

  // Setup algorithm chain: Input / Simulation / Output
  Simulation::setupInput(vars, sequencer, randomNumbers);
  setupGeant4Simulation(vars, sequencer, std::move(detectorConstructionFactory),
                        std::move(trackingGeometry), magneticField);
  Simulation::setupOutput(vars, sequencer);

  auto result = sequencer.run();
  return result;
}

}  // namespace ActsExamples
