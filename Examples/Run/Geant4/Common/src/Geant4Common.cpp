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
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"
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
#include <boost/program_options.hpp>

namespace ActsExamples {

void setupGeant4Simulation(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<G4RunManager> runManager,
    std::unique_ptr<G4VUserDetectorConstruction> detector,
    std::vector<G4UserRunAction*> runActions,
    std::vector<G4UserEventAction*> eventActions,
    std::vector<G4UserTrackingAction*> trackingActions,
    std::vector<G4UserSteppingAction*> steppingActions,
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
    const std::shared_ptr<const Acts::MagneticFieldProvider>& magneticField,
    bool materialRecording) {
  auto g4loglevel =
      Acts::Logging::Level(vars["g4-loglevel"].as<unsigned int>());

  size_t seed = vars["g4-seed"].as<size_t>();

  // Set the main Geant4 algorithm, primary generation, detector construction
  Geant4Simulation::Config g4Cfg;

  g4Cfg.runManager = std::move(runManager);
  g4Cfg.randomNumbers = std::make_shared<ActsExamples::RandomNumbers>(
      ActsExamples::RandomNumbers::Config{seed});
  g4Cfg.inputParticles = materialRecording ? Simulation::kParticlesInitial
                                           : Simulation::kParticlesSelection;

  // Read the particle from the generator
  SimParticleTranslation::Config g4PrCfg;

  if (materialRecording) {
    g4PrCfg.forcedPdgCode = 0;
    g4PrCfg.forcedCharge = 0.;
    g4PrCfg.forcedMass = 0.;
    // Set the material tracks at output
    g4Cfg.outputMaterialTracks = Simulation::kMaterialTracks;
  }

  // Set the primarty generator
  g4Cfg.primaryGeneratorAction = new SimParticleTranslation(
      g4PrCfg, Acts::getDefaultLogger("SimParticleTranslation", g4loglevel));
  g4Cfg.detectorConstruction = detector.release();

  // Set the user actions
  g4Cfg.runActions = std::move(runActions);
  g4Cfg.eventActions = std::move(eventActions);
  g4Cfg.trackingActions = std::move(trackingActions);
  g4Cfg.steppingActions = std::move(steppingActions);

  // An ACTS Magnetic field is provided
  if (magneticField) {
    MagneticFieldWrapper::Config g4FieldCfg;
    g4FieldCfg.magneticField = magneticField;
    g4Cfg.magneticField = new MagneticFieldWrapper(g4FieldCfg);
  }

  // An ACTS TrackingGeometry is provided, so simulation for sensitive
  // detectors is turned on - they need to get matched first
  if (trackingGeometry) {
    SensitiveSurfaceMapper::Config ssmCfg;
    ssmCfg.trackingGeometry = trackingGeometry;
    g4Cfg.sensitiveSurfaceMapper =
        std::make_shared<const SensitiveSurfaceMapper>(
            ssmCfg,
            Acts::getDefaultLogger("SensitiveSurfaceMapper", g4loglevel));

    // Add the output collection
    g4Cfg.outputSimHits = Simulation::kSimHits;
    g4Cfg.outputParticlesInitial = Simulation::kParticlesInitial;
    g4Cfg.outputParticlesFinal = Simulation::kParticlesFinal;
  }

  sequencer.addAlgorithm(std::make_shared<Geant4Simulation>(g4Cfg, g4loglevel));

  return;
}

int runMaterialRecording(
    const ActsExamples::Options::Variables& vars,
    std::unique_ptr<G4VUserDetectorConstruction> detector) {
  // Set up the sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));
  // Output level and log level
  auto logLevel = Options::readLogLevel(vars);
  auto g4loglevel =
      Acts::Logging::Level(vars["g4-loglevel"].as<unsigned int>());

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

  auto* physicsList = new MaterialPhysicsList(
      Acts::getDefaultLogger("MaterialPhysicsList", g4loglevel));

  auto runManager = std::make_shared<G4RunManager>();
  runManager->SetUserInitialization(physicsList);

  // The Geant4 actions needed
  std::vector<G4UserRunAction*> runActions = {};
  std::vector<G4UserEventAction*> eventActions = {};
  std::vector<G4UserTrackingAction*> trackingActions = {};

  MaterialSteppingAction::Config mStepCfg;
  mStepCfg.excludeMaterials = {"Air", "Vacuum"};
  std::vector<G4UserSteppingAction*> steppingActions = {
      new MaterialSteppingAction(
          mStepCfg,
          Acts::getDefaultLogger("MaterialSteppingAction", g4loglevel))};

  // Set up the Geant4 Simulation
  setupGeant4Simulation(vars, sequencer, runManager, std::move(detector),
                        runActions, eventActions, trackingActions,
                        steppingActions, nullptr, nullptr, true);

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
    std::unique_ptr<G4VUserDetectorConstruction> detector,
    const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry) {
  // Basic services
  auto randomNumbers =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));

  // Output level and log level
  auto g4loglevel =
      Acts::Logging::Level(vars["g4-loglevel"].as<unsigned int>());

  auto runManager = std::make_shared<G4RunManager>();
  runManager->SetUserInitialization(new FTFP_BERT());

  // The Geant4 actions needed
  std::vector<G4UserRunAction*> runActions = {};
  std::vector<G4UserEventAction*> eventActions = {};
  std::vector<G4UserTrackingAction*> trackingActions = {};
  std::vector<G4UserSteppingAction*> steppingActions = {};

  ParticleTrackingAction::Config g4TrackCfg;
  ParticleTrackingAction* particleAction = new ParticleTrackingAction(
      g4TrackCfg, Acts::getDefaultLogger("ParticleTrackingAction", g4loglevel));
  trackingActions.push_back(particleAction);

  SensitiveSteppingAction::Config g4StepCfg;
  SensitiveSteppingAction* sensitiveStepping = new SensitiveSteppingAction(
      g4StepCfg, Acts::getDefaultLogger("SensitiveSteppingAction", g4loglevel));
  steppingActions.push_back(sensitiveStepping);

  // Set up the sequencer
  Sequencer sequencer(Options::readSequencerConfig(vars));

  // Set up the magnetic field
  auto magneticField = ActsExamples::Options::readMagneticField(vars);

  // Setup algorithm chain: Input / Simulation / Output
  Simulation::setupInput(vars, sequencer, randomNumbers);
  setupGeant4Simulation(vars, sequencer, runManager, std::move(detector),
                        runActions, eventActions, trackingActions,
                        steppingActions, trackingGeometry, magneticField);
  Simulation::setupOutput(vars, sequencer);

  auto result = sequencer.run();
  return result;
}

}  // namespace ActsExamples
