// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Simulation.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/DetectorConstructionFactory.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/Geant4Manager.hpp"
#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"
#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"
#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"
#include "ActsExamples/Geant4/ParticleKillAction.hpp"
#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSteppingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"
#include "ActsExamples/Geant4/SimParticleTranslation.hpp"
#include "ActsExamples/Geant4/SteppingActionList.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>

#include <FTFP_BERT.hh>
#include <G4EmParameters.hh>
#include <G4FieldManager.hh>
#include <G4HadronicParameters.hh>
#include <G4HadronicProcessStore.hh>
#include <G4MagneticField.hh>
#include <G4RunManager.hh>
#include <G4StateManager.hh>
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

ActsExamples::Geant4SimulationBase::Geant4SimulationBase(
    const Config& cfg, std::string name, Acts::Logging::Level level)
    : IAlgorithm(std::move(name), level) {
  if (cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particle collection");
  }
  if (!cfg.detectorConstructionFactory) {
    throw std::invalid_argument("Missing detector construction factory");
  }
  if (!cfg.randomNumbers) {
    throw std::invalid_argument("Missing random numbers");
  }

  m_logger = Acts::getDefaultLogger("Geant4", level);

  m_eventStoreHolder = std::make_shared<EventStoreHolder>();

  // tweek logging
  {
    // If we are in VERBOSE mode, set the verbose level in Geant4 to 2.
    // 3 would be also possible, but that produces infinite amount of output.
    m_geant4Level = logger().level() == Acts::Logging::VERBOSE ? 2 : 0;
  }
}

ActsExamples::Geant4SimulationBase::~Geant4SimulationBase() = default;

void ActsExamples::Geant4SimulationBase::commonInitialization() {
  auto runManager = m_gean4Instance->runManager.get();

  // Set the detector construction
  {
    // Clear detector construction if it exists
    if (runManager->GetUserDetectorConstruction() != nullptr) {
      delete runManager->GetUserDetectorConstruction();
    }
    // G4RunManager will take care of deletion
    m_detectorConstruction =
        config().detectorConstructionFactory->factorize().release();
    runManager->SetUserInitialization(m_detectorConstruction);
    runManager->ReinitializeGeometry(true, true);
    runManager->InitializeGeometry();
  }
}

ActsExamples::ProcessCode ActsExamples::Geant4SimulationBase::initialize() {
  auto runManager = m_gean4Instance->runManager.get();

  // Initialize the Geant4 run manager
  runManager->Initialize();

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::Geant4SimulationBase::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Ensure exclusive access to the Geant4 run manager
  std::lock_guard<std::mutex> guard(m_gean4Instance->mutex);
  auto runManager = m_gean4Instance->runManager.get();

  // Set the seed new per event, so that we get reproducible results
  G4Random::setTheSeed(config().randomNumbers->generateSeed(ctx));

  // Get and reset event registry state
  auto& eventData = m_eventStoreHolder->store();
  eventData = EventStore{};

  // Register the current event store to the registry
  // this will allow access from the User*Actions
  eventData.store = &(ctx.eventStore);

  // Register the input particle read handle
  eventData.inputParticles = &m_inputParticles;

  ACTS_DEBUG("Sending Geant RunManager the BeamOn() command.");
  // Start simulation. each track is simulated as a separate Geant4 event.
  runManager->BeamOn(1);

  // Since these are std::set, this ensures that each particle is in both sets
  if (eventData.particlesInitial.size() != eventData.particlesFinal.size()) {
    ACTS_WARNING(
        "initial and final particle collections does not have the same size: "
        << eventData.particlesInitial.size() << " vs "
        << eventData.particlesFinal.size());
  }

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

  return ActsExamples::ProcessCode::SUCCESS;
}

std::shared_ptr<ActsExamples::Geant4Instance>
ActsExamples::Geant4SimulationBase::geant4Instance() const {
  return m_gean4Instance;
}

ActsExamples::Geant4Simulation::Geant4Simulation(const Config& cfg,
                                                 Acts::Logging::Level level)
    : Geant4SimulationBase(cfg, "Geant4Simulation", level), m_cfg(cfg) {
  m_gean4Instance =
      m_cfg.geant4Instance
          ? m_cfg.geant4Instance
          : Geant4Manager::instance().create(m_geant4Level, m_cfg.physicsList);
  if (m_gean4Instance->physicsListName != m_cfg.physicsList) {
    throw std::runtime_error("inconsistent physics list");
  }
  auto runManager = m_gean4Instance->runManager.get();

  commonInitialization();

  // Set the primarty generator
  {
    SimParticleTranslation::Config prCfg;
    prCfg.eventStoreHolder = m_eventStoreHolder;
    // G4RunManager will take care of deletion
    auto primaryGeneratorAction = new SimParticleTranslation(
        prCfg, m_logger->cloneWithSuffix("SimParticleTranslation"));
    // Clear primary generation action if it exists
    if (runManager->GetUserPrimaryGeneratorAction() != nullptr) {
      delete runManager->GetUserPrimaryGeneratorAction();
    }
    // Set the primary generator action
    runManager->SetUserAction(primaryGeneratorAction);
  }

  // Please note:
  // The following two blocks rely on the fact that the Acts
  // detector constructions cache the world volume

  // Get the g4World cache
  G4VPhysicalVolume* g4World = m_detectorConstruction->Construct();

  // Particle action
  {
    ParticleTrackingAction::Config trackingCfg;
    trackingCfg.eventStoreHolder = m_eventStoreHolder;
    trackingCfg.keepParticlesWithoutHits = cfg.keepParticlesWithoutHits;
    // G4RunManager will take care of deletion
    auto trackingAction = new ParticleTrackingAction(
        trackingCfg, m_logger->cloneWithSuffix("ParticleTracking"));
    // Clear tracking action if it exists
    if (runManager->GetUserTrackingAction() != nullptr) {
      delete runManager->GetUserTrackingAction();
    }
    runManager->SetUserAction(trackingAction);
  }

  // Stepping actions
  {
    SensitiveSteppingAction::Config stepCfg;
    stepCfg.eventStoreHolder = m_eventStoreHolder;
    stepCfg.charged = true;
    stepCfg.neutral = false;
    stepCfg.primary = true;
    stepCfg.secondary = cfg.recordHitsOfSecondaries;

    ParticleKillAction::Config particleKillCfg;
    particleKillCfg.volume = cfg.killVolume;
    particleKillCfg.maxTime = cfg.killAfterTime;

    SteppingActionList::Config steppingCfg;
    steppingCfg.actions.push_back(std::make_unique<SensitiveSteppingAction>(
        stepCfg, m_logger->cloneWithSuffix("SensitiveStepping")));
    steppingCfg.actions.push_back(std::make_unique<ParticleKillAction>(
        particleKillCfg, m_logger->cloneWithSuffix("Killer")));
    // G4RunManager will take care of deletion
    auto steppingAction = new SteppingActionList(steppingCfg);
    // Clear stepping action if it exists
    if (runManager->GetUserSteppingAction() != nullptr) {
      delete runManager->GetUserSteppingAction();
    }
    runManager->SetUserAction(steppingAction);
  }

  // Set the magnetic field
  if (cfg.magneticField) {
    ACTS_INFO("Setting ACTS configured field to Geant4.");

    MagneticFieldWrapper::Config g4FieldCfg;
    g4FieldCfg.magneticField = cfg.magneticField;
    m_magneticField = std::make_unique<MagneticFieldWrapper>(g4FieldCfg);

    // Set the field ot the G4Field manager
    m_fieldManager = std::make_unique<G4FieldManager>();
    m_fieldManager->SetDetectorField(m_magneticField.get());
    m_fieldManager->CreateChordFinder(m_magneticField.get());

    // Propagate down to all childrend
    g4World->GetLogicalVolume()->SetFieldManager(m_fieldManager.get(), true);
  }

  // An ACTS TrackingGeometry is provided, so simulation for sensitive
  // detectors is turned on - they need to get matched first
  if (cfg.trackingGeometry) {
    ACTS_INFO(
        "Remapping selected volumes from Geant4 to Acts::Surface::GeometryID");

    SensitiveSurfaceMapper::Config ssmCfg;
    ssmCfg.trackingGeometry = cfg.trackingGeometry;
    ssmCfg.volumeMappings = cfg.volumeMappings;
    ssmCfg.materialMappings = cfg.materialMappings;

    SensitiveSurfaceMapper sensitiveSurfaceMapper(
        ssmCfg, m_logger->cloneWithSuffix("SensitiveSurfaceMapper"));
    int sCounter = 0;
    sensitiveSurfaceMapper.remapSensitiveNames(
        g4World, Acts::Transform3::Identity(), sCounter);

    ACTS_INFO("Remapping successful for " << sCounter << " selected volumes.");
  }

  m_inputParticles.initialize(cfg.inputParticles);
  m_outputSimHits.initialize(cfg.outputSimHits);
  m_outputParticlesInitial.initialize(cfg.outputParticlesInitial);
  m_outputParticlesFinal.initialize(cfg.outputParticlesFinal);
}

ActsExamples::Geant4Simulation::~Geant4Simulation() = default;

ActsExamples::ProcessCode ActsExamples::Geant4Simulation::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  Geant4SimulationBase::execute(ctx);

  auto& eventData = m_eventStoreHolder->store();

  // Output handling: Simulation
  m_outputParticlesInitial(
      ctx, SimParticleContainer(eventData.particlesInitial.begin(),
                                eventData.particlesInitial.end()));
  m_outputParticlesFinal(ctx,
                         SimParticleContainer(eventData.particlesFinal.begin(),
                                              eventData.particlesFinal.end()));

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

  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::Geant4MaterialRecording::Geant4MaterialRecording(
    const Config& cfg, Acts::Logging::Level level)
    : Geant4SimulationBase(cfg, "Geant4Simulation", level), m_cfg(cfg) {
  auto physicsListName = "MaterialPhysicsList";
  m_gean4Instance =
      m_cfg.geant4Instance
          ? m_cfg.geant4Instance
          : Geant4Manager::instance().create(
                m_geant4Level,
                std::make_unique<MaterialPhysicsList>(
                    m_logger->cloneWithSuffix("MaterialPhysicsList")),
                physicsListName);
  if (m_gean4Instance->physicsListName != physicsListName) {
    throw std::runtime_error("inconsistent physics list");
  }
  auto runManager = m_gean4Instance->runManager.get();

  commonInitialization();

  // Set the primarty generator
  {
    SimParticleTranslation::Config prCfg;
    prCfg.eventStoreHolder = m_eventStoreHolder;
    prCfg.forcedPdgCode = 0;
    prCfg.forcedCharge = 0.;
    prCfg.forcedMass = 0.;
    // G4RunManager will take care of deletion
    auto primaryGeneratorAction = new SimParticleTranslation(
        prCfg, m_logger->cloneWithSuffix("SimParticleTranslation"));
    // Clear primary generation action if it exists
    if (runManager->GetUserPrimaryGeneratorAction() != nullptr) {
      delete runManager->GetUserPrimaryGeneratorAction();
    }
    // Set the primary generator action
    runManager->SetUserAction(primaryGeneratorAction);
  }

  m_inputParticles.initialize(cfg.inputParticles);
  m_outputMaterialTracks.initialize(cfg.outputMaterialTracks);
}

ActsExamples::Geant4MaterialRecording::~Geant4MaterialRecording() = default;

ActsExamples::ProcessCode ActsExamples::Geant4MaterialRecording::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  Geant4SimulationBase::execute(ctx);

  auto& eventData = m_eventStoreHolder->store();

  // Output handling: Material tracks
  m_outputMaterialTracks(
      ctx, decltype(eventData.materialTracks)(eventData.materialTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}
