// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/InteractionProcessRecording.hpp"
#include <iostream>
#include <stdexcept>
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "OREventAction.hpp"
#include "ORPrimaryGeneratorAction.hpp"
#include "ORSteppingAction.hpp"
#include "RunAction.hpp"
#include "FTFP_BERT.hh"
#include "ActsExamples/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"


ActsExamples::InteractionProcessRecording::InteractionProcessRecording(
    const ActsExamples::InteractionProcessRecording::Config& cnf,
    Acts::Logging::Level                 level)
  : ActsExamples::BareAlgorithm("InteractionProcessRecording", level)
  , m_cfg(cnf)
  , m_runManager(std::make_unique<G4RunManager>())
{
  /// Check if the geometry should be accessed over the geant4 service
  if (m_cfg.geant4Service) {
    m_runManager->SetUserInitialization(m_cfg.geant4Service->geant4Geometry());
  } else if (!m_cfg.gdmlFile.empty()) {
    /// Access the geometry from the gdml file
    ACTS_INFO(
        "received Geant4 geometry from GDML file: " << m_cfg.gdmlFile.c_str());
    ActsExamples::GdmlDetectorConstruction* detConstruction
        = new ActsExamples::GdmlDetectorConstruction();
    detConstruction->setGdmlInput(m_cfg.gdmlFile.c_str());
    m_runManager->SetUserInitialization(
        detConstruction);  // constructs detector (calls Construct in
                           // Geant4DetectorConstruction)
  } else {
    throw std::invalid_argument("Missing geometry input for Geant4");
  }

  /// Now set up the Geant4 simulation
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new ActsExamples::Geant4::ORPrimaryGeneratorAction(
      cnf.pdg, cnf.momentum * 1000., cnf.lockAngle, cnf.phi, cnf.theta, cnf.lockPosition, {cnf.pos.x(), cnf.pos.y(), cnf.pos.z()}, m_cfg.seed1, m_cfg.seed2));
  ActsExamples::RunAction* runaction = new ActsExamples::MMRunAction();
  m_runManager->SetUserAction(runaction);
  ActsExamples::OREventAction* evtAct = new ActsExamples::OREventAction();
  m_runManager->SetUserAction(evtAct);
  m_runManager->SetUserAction(new ActsExamples::ORSteppingAction(evtAct));
  m_runManager->Initialize();
}

ActsExamples::ProcessCode
ActsExamples::InteractionProcessRecording::execute(const ActsExamples::AlgorithmContext& context) const
{
  // Begin with the simulation
  m_runManager->BeamOn(m_cfg.tracksPerEvent);
  // Retrieve the track material tracks from Geant4
  auto recordedParticles
      = ActsExamples::OREventAction::Instance()->outcomingParticles(m_cfg.pdg, m_cfg.momentum, m_cfg.phi, m_cfg.theta); // Keeping momentum in Acts units
  ACTS_INFO("Received " << recordedParticles.particles.size()
                        << " particles. Writing them now onto file...");

  // Write the recorded material to the event store
  context.eventStore.add(m_cfg.particleCollection,
                         std::move(recordedParticles));

  return ActsExamples::ProcessCode::SUCCESS;
}
