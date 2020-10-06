// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "OREventAction.hpp"
#include "ORPrimaryGeneratorAction.hpp"
#include "ORSteppingAction.hpp"
#include "RunAction.hpp"
#include "ActsExamples/Geant4/InteractionProcessRecording.hpp"
#include <iostream>
#include <stdexcept>
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
//~ #include "OREventAction.hpp"
//~ #include "ORPrimaryGeneratorAction.hpp"
//~ #include "ORSteppingAction.hpp"
//~ #include "RunAction.hpp"
#include "FTFP_BERT.hh"
#include "ActsExamples/EventData/SimParticle.hpp"

ActsExamples::InteractionProcessRecording::InteractionProcessRecording(
    ActsExamples::InteractionProcessRecording::Config&& cnf,
    Acts::Logging::Level                 level)
  : ActsExamples::BareAlgorithm("InteractionProcessRecording", level)
  , m_cfg(std::move(cnf))
  , m_runManager(std::make_unique<G4RunManager>())
{
	if(m_cfg.eventInput.empty())
	{
		throw std::invalid_argument("Missing input event collection");
	}
	if(m_cfg.eventOutput.empty())
	{
		throw std::invalid_argument("Missing output event collection");
	}
	if(!m_cfg.detectorConstruction)
	{
		throw std::invalid_argument("Missing detector construction object");
	}

  /// Now set up the Geant4 simulation
  m_runManager->SetUserInitialization(m_cfg.detectorConstruction.release());
  m_runManager->SetUserInitialization(new FTFP_BERT);
  //~ m_runManager->SetUserAction(new ActsExamples::ORPrimaryGeneratorAction(
      //~ cnf.pdg, cnf.momentum * 1000., cnf.lockAngle, cnf.phi, cnf.theta, cnf.lockPosition, {cnf.pos.x(), cnf.pos.y(), cnf.pos.z()}, m_cfg.seed1, m_cfg.seed2));
  m_runManager->SetUserAction(new ActsExamples::ORPrimaryGeneratorAction(m_cfg.seed1, m_cfg.seed2));
  m_runManager->SetUserAction(new ActsExamples::RunAction());
  m_runManager->SetUserAction(new ActsExamples::OREventAction());
  m_runManager->SetUserAction(new ActsExamples::ORSteppingAction());
  m_runManager->Initialize();
}

ActsExamples::ProcessCode
ActsExamples::InteractionProcessRecording::execute(const ActsExamples::AlgorithmContext& context) const
{
  // ensure exclusive access to the geant run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);
  
  // TODO: params lesen
  const auto initialParticles = context.eventStore.get<ActsExamples::SimParticleContainer>(m_cfg.eventInput);
std::cout << typeid(initialParticles).name() << std::endl;
      
  //~ ActsExamples::ORPrimaryGeneratorAction::instance()->prepareParticleGun(G4int pdg,
     //~ G4double        momentum,
	 //~ G4double phi,
	 //~ G4double theta,
	 //~ G4ThreeVector pos);
  
  // Begin with the simulation
  m_runManager->BeamOn(m_cfg.tracksPerEvent);
    
  // Retrieve the track material tracks from Geant4
  //~ auto recordedParticles
      //~ = ActsExamples::OREventAction::instance()->processTracks(); // TODO: might need to store initial parameters additionally | m_cfg.pdg, m_cfg.momentum, m_cfg.phi, m_cfg.theta); // Keeping momentum in Acts units
  
  
  //~ ACTS_INFO("Received " << recordedParticles.particles.size()
                        //~ << " particles. Writing them now onto file...");

  //~ // Write the recorded material to the event store
  //~ context.eventStore.add(m_cfg.eventOutput,
                         //~ std::move(recordedParticles));

  return ActsExamples::ProcessCode::SUCCESS;
}
