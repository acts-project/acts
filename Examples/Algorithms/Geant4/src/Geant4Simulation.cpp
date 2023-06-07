// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Simulation.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/ActsSteppingActionList.hpp"
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"
#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"
#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"
#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"
#include "ActsExamples/Geant4/ParticleKillAction.hpp"
#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSteppingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include <iostream>
#include <stdexcept>

#include <FTFP_BERT.hh>
#include <G4EmParameters.hh>
#include <G4FieldManager.hh>
#include <G4HadronicParameters.hh>
#include <G4HadronicProcessStore.hh>
#include <G4MagneticField.hh>
#include <G4RunManager.hh>
#include <G4TransportationManager.hh>
#include <G4UniformMagField.hh>
#include <G4UserEventAction.hh>
#include <G4UserLimits.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4Version.hh>

ActsExamples::Geant4Simulation::Config
ActsExamples::Geant4Simulation::makeGeant4Config(
    const Acts::Logger& logger,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    std::shared_ptr<G4VUserDetectorConstruction> detector,
    std::shared_ptr<G4VUserPhysicsList> physicsList,
    const SimParticleTranslation::Config& prCfg) {
  Geant4Simulation::Config g4Cfg;

  // Set the main Geant4 algorithm, primary generation, detector
  // construction
  g4Cfg.randomNumbers = std::move(randomNumbers);
  g4Cfg.runManager = std::make_shared<G4RunManager>();
  g4Cfg.runManager->SetUserInitialization(physicsList.get());

  // Set the primarty generator
  g4Cfg.primaryGeneratorAction = std::make_shared<SimParticleTranslation>(
      prCfg, logger.cloneWithSuffix("SimParticleTranslation"));
  g4Cfg.detectorConstruction = detector;

  return g4Cfg;
}

ActsExamples::Geant4Simulation::Config
ActsExamples::Geant4Simulation::makeGeant4MaterialRecordingConfig(
    Acts::Logging::Level level,
    std::shared_ptr<G4VUserDetectorConstruction> detector,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    const std::string& inputParticles,
    const std::string& outputMaterialTracks) {
  auto logger = Acts::getDefaultLogger("Geant4", level);
  auto physicsList = std::make_shared<MaterialPhysicsList>(
      logger->cloneWithSuffix("MaterialPhysicsList"));

  // Read the particle from the generator
  SimParticleTranslation::Config g4PrCfg;
  g4PrCfg.forcedPdgCode = 0;
  g4PrCfg.forcedCharge = 0.;
  g4PrCfg.forcedMass = 0.;

  auto g4Cfg = makeGeant4Config(*logger, std::move(randomNumbers),
                                std::move(detector), physicsList, g4PrCfg);
  g4Cfg.inputParticles = inputParticles;

  MaterialSteppingAction::Config mStepCfg;
  mStepCfg.excludeMaterials = {"Air", "Vacuum"};
  auto steppingAction = std::make_shared<MaterialSteppingAction>(
      mStepCfg, logger->cloneWithSuffix("MaterialSteppingAction"));
  g4Cfg.steppingAction = steppingAction;

  // Set the material tracks at output
  g4Cfg.outputMaterialTracks = outputMaterialTracks;

  return g4Cfg;
}

ActsExamples::Geant4Simulation::Config
ActsExamples::Geant4Simulation::makeGeant4SimulationConfig(
    Acts::Logging::Level level,
    std::shared_ptr<G4VUserDetectorConstruction> detector,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    const std::string& inputParticles,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    const std::vector<std::string>& volumeMappings,
    const std::vector<std::string>& materialMappings,
    std::shared_ptr<const Acts::Volume> killVolume, double killAfterTime,
    bool recordHitsOfSecondaries, bool keepParticlesWithoutHits) {
  auto logger = Acts::getDefaultLogger("Geant4", level);

  auto physicsList = std::make_shared<FTFP_BERT>();
  auto g4Cfg = makeGeant4Config(*logger, std::move(randomNumbers),
                                std::move(detector), std::move(physicsList),
                                SimParticleTranslation::Config{});
  g4Cfg.inputParticles = inputParticles;

  // Particle action
  ParticleTrackingAction::Config trackingCfg;
  trackingCfg.keepParticlesWithoutHits = keepParticlesWithoutHits;
  g4Cfg.trackingAction = std::make_shared<ParticleTrackingAction>(
      trackingCfg, logger->cloneWithSuffix("ParticleTracking"));

  // Stepping actions
  ActsSteppingActionList::Config steppingCfg;

  SensitiveSteppingAction::Config g4StepCfg;
  g4StepCfg.charged = true;
  g4StepCfg.neutral = false;
  g4StepCfg.primary = true;
  g4StepCfg.secondary = recordHitsOfSecondaries;
  steppingCfg.actions.push_back(new SensitiveSteppingAction(
      g4StepCfg, logger->cloneWithSuffix("SensitiveStepping")));

  steppingCfg.actions.push_back(new ParticleKillAction(
      ParticleKillAction::Config{killVolume, killAfterTime},
      logger->cloneWithSuffix("Killer")));

  g4Cfg.steppingAction = std::make_shared<ActsSteppingActionList>(steppingCfg);

  // An ACTS Magnetic field is provided
  if (magneticField) {
    MagneticFieldWrapper::Config g4FieldCfg;
    g4FieldCfg.magneticField = magneticField;
    g4Cfg.magneticField = std::make_shared<MagneticFieldWrapper>(g4FieldCfg);
  }

  // An ACTS TrackingGeometry is provided, so simulation for sensitive
  // detectors is turned on - they need to get matched first
  if (trackingGeometry) {
    SensitiveSurfaceMapper::Config ssmCfg;
    ssmCfg.trackingGeometry = trackingGeometry;

    // Take the default args if nothing provided
    if (not volumeMappings.empty()) {
      ssmCfg.volumeMappings = volumeMappings;
    }
    if (not materialMappings.empty()) {
      ssmCfg.materialMappings = materialMappings;
    }

    g4Cfg.sensitiveSurfaceMapper =
        std::make_shared<const SensitiveSurfaceMapper>(
            ssmCfg, logger->cloneWithSuffix("SensitiveSurfaceMapper"));
  }

  return g4Cfg;
}

ActsExamples::Geant4Simulation::Geant4Simulation(
    const ActsExamples::Geant4Simulation::Config& config,
    Acts::Logging::Level level)
    : IAlgorithm("Geant4Simulation", level), m_cfg(config) {
  if (m_cfg.detectorConstruction == nullptr) {
    throw std::invalid_argument("Missing G4 DetectorConstruction object");
  }
  if (m_cfg.primaryGeneratorAction == nullptr) {
    throw std::invalid_argument("Missing G4 PrimaryGeneratorAction object");
  }
  if (!m_cfg.runManager) {
    throw std::invalid_argument("Missing G4 RunManager object");
  }
  if (!m_cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers tool");
  }

  if (m_cfg.sensitiveSurfaceMapper) {
    if (m_cfg.outputSimHits.empty()) {
      ACTS_WARNING("No output sim hits collection configured");
    }
    m_outputSimHits.initialize(m_cfg.outputSimHits);

    if (m_cfg.outputParticlesInitial.empty()) {
      ACTS_WARNING("No output initial particles collection configured");
    }
    m_outputParticlesInitial.initialize(m_cfg.outputParticlesInitial);

    if (m_cfg.outputParticlesFinal.empty()) {
      ACTS_WARNING("No output final particles collection configured");
    }
    m_outputParticlesFinal.initialize(m_cfg.outputParticlesFinal);
  }

  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }

  // If we are in VERBOSE mode, set the verbose level in Geant4 to 2.
  // 3 would be also possible, but that produces infinite amount of output.
  const int geantVerboseLevel =
      logger().level() == Acts::Logging::VERBOSE ? 2 : 0;
  m_cfg.runManager->SetVerboseLevel(geantVerboseLevel);
  G4EventManager::GetEventManager()->SetVerboseLevel(geantVerboseLevel);
  G4EventManager::GetEventManager()->GetTrackingManager()->SetVerboseLevel(
      geantVerboseLevel);
  G4EventManager::GetEventManager()->GetStackManager()->SetVerboseLevel(
      geantVerboseLevel);

  // Suppress the printing of physics information.
#if G4VERSION_NUMBER >= 1100
  G4HadronicParameters::Instance()->SetVerboseLevel(0);
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  G4EmParameters::Instance()->SetIsPrintedFlag(true);
#endif

  // Set the detector construction
  m_cfg.runManager->SetUserInitialization(m_cfg.detectorConstruction.get());

  // Set the primary generator action
  m_cfg.runManager->SetUserAction(m_cfg.primaryGeneratorAction.get());

  // Set the configured user actions
  m_cfg.runManager->SetUserAction(m_cfg.runAction.get());
  m_cfg.runManager->SetUserAction(m_cfg.eventAction.get());
  m_cfg.runManager->SetUserAction(m_cfg.trackingAction.get());
  m_cfg.runManager->SetUserAction(m_cfg.steppingAction.get());

  // Initialize the Geant4 run manager
  m_cfg.runManager->Initialize();

  // Please note:
  // The following two blocks rely on the fact that the Acts
  // detector constructions cache the world volume

  // Set the magnetic field
  if (m_cfg.magneticField != nullptr) {
    ACTS_INFO("Setting ACTS configured field to Geant4.");
    // Get the g4World cache
    G4VPhysicalVolume* g4World = m_cfg.detectorConstruction->Construct();
    /// Set the field ot the G4Field manager
    G4FieldManager* fieldMgr = new G4FieldManager();
    fieldMgr->SetDetectorField(m_cfg.magneticField.get());
    fieldMgr->CreateChordFinder(m_cfg.magneticField.get());

    // Propagate down to all childrend
    g4World->GetLogicalVolume()->SetFieldManager(fieldMgr, true);
  }

  // Map simulation to reconstruction geometry
  if (m_cfg.sensitiveSurfaceMapper != nullptr) {
    ACTS_INFO(
        "Remapping selected volumes from Geant4 to Acts::Surface::GeometryID");

    G4VPhysicalVolume* g4World = m_cfg.detectorConstruction->Construct();
    int sCounter = 0;
    m_cfg.sensitiveSurfaceMapper->remapSensitiveNames(
        g4World, Acts::Transform3::Identity(), sCounter);

    ACTS_INFO("Remapping successful for " << sCounter << " selected volumes.");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputMaterialTracks.maybeInitialize(m_cfg.outputMaterialTracks);
}

ActsExamples::Geant4Simulation::~Geant4Simulation() = default;

ActsExamples::ProcessCode ActsExamples::Geant4Simulation::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Ensure exclusive access to the Geant4 run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  // Set the seed new per event, so that we get reproducible results
  G4Random::setTheSeed(m_cfg.randomNumbers->generateSeed(ctx));

  // Get and reset event registry state
  auto& eventData = EventStoreRegistry::eventData();
  eventData = EventStoreRegistry::State{};

  // Register the current event store to the registry
  // this will allow access from the User*Actions
  eventData.store = &(ctx.eventStore);

  // Register the input particle read handle
  eventData.inputParticles = &m_inputParticles;

  ACTS_DEBUG("Sending Geant RunManager the BeamOn() command.");
  // Start simulation. each track is simulated as a separate Geant4 event.
  m_cfg.runManager->BeamOn(1);

  // Since these are std::set, this ensures that each particle is in both sets
  assert(eventData.particlesInitial.size() == eventData.particlesFinal.size());

  // Print out warnings about possible particle collision if happened
  if (eventData.particleIdCollisionsInitial > 0 or
      eventData.particleIdCollisionsFinal > 0 or
      eventData.parentIdNotFound > 0) {
    ACTS_WARNING(
        "Particle ID collisions detected, don't trust the particle "
        "identification!");
    ACTS_WARNING("- initial states: " << eventData.particleIdCollisionsInitial);
    ACTS_WARNING("- final states: " << eventData.particleIdCollisionsFinal);
    ACTS_WARNING("- parent ID not found: " << eventData.parentIdNotFound);
  }

  // Output handling: Initial/Final particles
  if (not m_cfg.outputParticlesInitial.empty() and
      not m_cfg.outputParticlesFinal.empty()) {
    m_outputParticlesInitial(
        ctx, SimParticleContainer(eventData.particlesInitial.begin(),
                                  eventData.particlesInitial.end()));
    m_outputParticlesFinal(
        ctx, SimParticleContainer(eventData.particlesFinal.begin(),
                                  eventData.particlesFinal.end()));
  }

  // Output handling: Simulated hits
  if (not m_cfg.outputSimHits.empty()) {
#if BOOST_VERSION < 107800
    SimHitContainer container;
    for (const auto& hit : eventData.hits) {
      container.insert(hit);
    }
    m_outputSimHits(ctx, std::move(container));
#else
    m_outputSimHits(
        ctx, SimHitContainer(eventData.hits.begin(), eventData.hits.end()));
#endif
  }

  // Output handling: Material tracks
  if (not m_cfg.outputMaterialTracks.empty()) {
    m_outputMaterialTracks(
        ctx, decltype(eventData.materialTracks)(eventData.materialTracks));
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
