// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/GeantinoRecording.hpp"

#include <FTFP_BERT.hh>
#include <G4RunManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <iostream>
#include <stdexcept>

#include "ACTFW/Framework/WhiteBoard.hpp"
#include "EventAction.hpp"
#include "PrimaryGeneratorAction.hpp"
#include "RunAction.hpp"
#include "SteppingAction.hpp"

using namespace ActsExamples;

GeantinoRecording::GeantinoRecording(GeantinoRecording::Config&& cfg,
                                     Acts::Logging::Level lvl)
    : BareAlgorithm("GeantinoRecording", lvl),
      m_cfg(std::move(cfg)),
      m_runManager(std::make_unique<G4RunManager>()) {
  if (m_cfg.outputMaterialTracks.empty()) {
    throw std::invalid_argument("Missing output material tracks collection");
  }
  if (not m_cfg.detectorConstruction) {
    throw std::invalid_argument("Missing detector construction object");
  }

  m_runManager->SetUserInitialization(m_cfg.detectorConstruction.release());
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new RunAction());
  m_runManager->SetUserAction(new EventAction());
  m_runManager->SetUserAction(
      new PrimaryGeneratorAction("geantino", 1000., m_cfg.seed1, m_cfg.seed2));
  m_runManager->SetUserAction(new SteppingAction());
  m_runManager->Initialize();
}

// needed to allow std::unique_ptr<G4RunManager> with forward-declared class.
GeantinoRecording::~GeantinoRecording() {}

FW::ProcessCode GeantinoRecording::execute(
    const FW::AlgorithmContext& ctx) const {
  // ensure exclusive access to the geant run manager
  std::lock_guard<std::mutex> guard(m_runManagerLock);

  // TODO use framework random numbers directly or at least context seed
  // TODO take particles collection as input instead of generating them

  // start simulation. each track is simulated as a separate Geant4 event.
  m_runManager->BeamOn(m_cfg.tracksPerEvent);

  auto materialTracks = EventAction::instance()->materialTracks();
  // Write the recorded material to the event store
  ctx.eventStore.add(m_cfg.outputMaterialTracks, move(materialTracks));

  return FW::ProcessCode::SUCCESS;
}
