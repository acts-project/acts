// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Plugins/HepMC3/EventRecording.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include <iostream>
#include <stdexcept>
#include "EventAction.hpp"
#include "FTFP_BERT.hh"
#include "PrimaryGeneratorAction.hpp"
#include "RunAction.hpp"
#include "SteppingAction.hpp"

#include <HepMC3/GenEvent.h>

ActsExamples::EventRecording::EventRecording(
    ActsExamples::EventRecording::Config&& cnf, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("EventRecording", level),
      m_cfg(std::move(cnf)),
      m_runManager(std::make_unique<G4RunManager>()) {
  if (m_cfg.eventInput.empty()) {
    throw std::invalid_argument("Missing input event collection");
  }
  if (m_cfg.eventOutput.empty()) {
    throw std::invalid_argument("Missing output event collection");
  }
  if (!m_cfg.detectorConstruction) {
    throw std::invalid_argument("Missing detector construction object");
  }

  /// Now set up the Geant4 simulation
  m_runManager->SetUserInitialization(m_cfg.detectorConstruction.release());
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new ActsExamples::RunAction());
  m_runManager->SetUserAction(new ActsExamples::EventAction());
  m_runManager->SetUserAction(
      new ActsExamples::PrimaryGeneratorAction(m_cfg.seed1, m_cfg.seed2));
  m_runManager->SetUserAction(new ActsExamples::SteppingAction());
  m_runManager->Initialize();
}

ActsExamples::ProcessCode ActsExamples::EventRecording::execute(
    const ActsExamples::AlgorithmContext& context) const {
  // ensure exclusive access to the geant run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  // Retrieve the initial particles
  const auto initialParticles =
      context.eventStore.get<ActsExamples::SimParticleContainer>(
          m_cfg.eventInput);

  // Storage of events that will be produced
  std::vector<std::shared_ptr<HepMC3::GenEvent>> events;
  events.reserve(initialParticles.size());

  for (const auto& part : initialParticles) {
    // Prepare the particle gun
    const auto pos = part.position();
    const auto dir = part.unitDirection();
    ActsExamples::PrimaryGeneratorAction::instance()->prepareParticleGun(
        part.pdg(), part.absMomentum(), {pos[0], pos[1], pos[2]},
        {dir[0], dir[1], dir[2]});

    // Begin with the simulation
    m_runManager->BeamOn(1);

    // Store the result
    events.push_back(ActsExamples::EventAction::instance()->event());
  }

  ACTS_INFO(events.size() << " events generated");

  // Write the recorded material to the event store
  context.eventStore.add(m_cfg.eventOutput, std::move(events));

  return ActsExamples::ProcessCode::SUCCESS;
}
