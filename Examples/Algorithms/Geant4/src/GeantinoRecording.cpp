// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/GeantinoRecording.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"

#include <iostream>
#include <stdexcept>

#include <FTFP_BERT.hh>
#include <G4RunManager.hh>
#include <G4VUserDetectorConstruction.hh>

#include "EventAction.hpp"
#include "RunAction.hpp"
#include "SteppingAction.hpp"

using namespace ActsExamples;
using namespace ActsExamples::Geant4;

GeantinoRecording::GeantinoRecording(GeantinoRecording::Config config,
                                     Acts::Logging::Level level)
    : BareAlgorithm("GeantinoRecording", level),
      m_cfg(std::move(config)),
      m_runManager(std::make_unique<G4RunManager>()) {
  std::cout << __FILE__ << ":" << __LINE__ << std::endl;
  if (m_cfg.outputMaterialTracks.empty()) {
    throw std::invalid_argument("Missing output material tracks collection");
  }
  m_cfg.generationConfig.particleName = "geantino";
  m_cfg.generationConfig.energy = 1000.;

  if (!m_cfg.detectorConstructionFactory) {
    throw std::invalid_argument("Missing detector construction object factory");
  }

  // This object here retains owner
  m_runManager->SetUserInitialization(
      (*m_cfg.detectorConstructionFactory)().release());
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new RunAction());
  m_runManager->SetUserAction(new EventAction());
  m_runManager->SetUserAction(
      new PrimaryGeneratorAction(m_cfg.generationConfig));
  m_runManager->SetUserAction(new SteppingAction());
  m_runManager->Initialize();
}

// needed to allow std::unique_ptr<G4RunManager> with forward-declared class.
GeantinoRecording::~GeantinoRecording() {}

ActsExamples::ProcessCode GeantinoRecording::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // ensure exclusive access to the geant run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  // TODO use framework random numbers directly or at least context seed
  // TODO take particles collection as input instead of generating them

  // start simulation. each track is simulated as a separate Geant4 event.
  m_runManager->BeamOn(m_cfg.tracksPerEvent);

  auto materialTracks = EventAction::instance()->materialTracks();
  // Write the recorded material to the event store
  ctx.eventStore.add(m_cfg.outputMaterialTracks, move(materialTracks));

  return ActsExamples::ProcessCode::SUCCESS;
}
