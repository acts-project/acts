// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/GeantinoRecording.hpp"

#include <FTFP_BERT.hh>
#include <G4RunManager.hh>
#include <iostream>
#include <stdexcept>

#include "ACTFW/Framework/WhiteBoard.hpp"
#include "MMDetectorConstruction.hpp"
#include "MMEventAction.hpp"
#include "MMPrimaryGeneratorAction.hpp"
#include "MMRunAction.hpp"
#include "MMSteppingAction.hpp"

FW::GeantinoRecording::GeantinoRecording(
    const FW::GeantinoRecording::Config& cfg, Acts::Logging::Level lvl)
    : FW::BareAlgorithm("GeantinoRecording", lvl),
      m_cfg(cfg),
      m_runManager(std::make_unique<G4RunManager>()) {
  /// Check if the geometry should be accessed over the geant4 service
  if (m_cfg.geant4Service) {
    m_runManager->SetUserInitialization(m_cfg.geant4Service->geant4Geometry());
  } else if (!m_cfg.gdmlFile.empty()) {
    /// Access the geometry from the gdml file
    ACTS_INFO(
        "received Geant4 geometry from GDML file: " << m_cfg.gdmlFile.c_str());
    FW::Geant4::MMDetectorConstruction* detConstruction =
        new FW::Geant4::MMDetectorConstruction();
    detConstruction->setGdmlInput(m_cfg.gdmlFile.c_str());
    m_runManager->SetUserInitialization(
        detConstruction);  // constructs detector (calls Construct in
                           // Geant4DetectorConstruction)
  } else {
    throw std::invalid_argument("Missing geometry input for Geant4");
  }

  /// Now set up the Geant4 simulation
  m_runManager->SetUserInitialization(new FTFP_BERT);
  m_runManager->SetUserAction(new FW::Geant4::MMPrimaryGeneratorAction(
      "geantino", 1000., m_cfg.seed1, m_cfg.seed2));
  FW::Geant4::MMRunAction* runaction = new FW::Geant4::MMRunAction();
  m_runManager->SetUserAction(runaction);
  m_runManager->SetUserAction(new FW::Geant4::MMEventAction());
  m_runManager->SetUserAction(new FW::Geant4::MMSteppingAction());
  m_runManager->Initialize();
}

// needed to allow std::unique_ptr<G4RunManager> with forward-declared class.
FW::GeantinoRecording::~GeantinoRecording() {}

FW::ProcessCode FW::GeantinoRecording::execute(
    const FW::AlgorithmContext& context) const {
  // Begin with the simulation
  m_runManager->BeamOn(m_cfg.tracksPerEvent);
  // Retrieve the track material tracks from Geant4
  auto recordedMaterial =
      FW::Geant4::MMEventAction::Instance()->MaterialTracks();
  ACTS_INFO("Received " << recordedMaterial.size()
                        << " MaterialTracks. Writing them now onto file...");

  // Write the recorded material to the event store
  context.eventStore.add(m_cfg.geantMaterialCollection,
                         std::move(recordedMaterial));

  // // Retrieve the sim hit track steps from Geant4
  // auto trackSteps = FW::Geant4::MMEventAction::Instance()->TrackSteps();
  // ACTS_INFO("Received " << trackSteps.size()
  //                       << " steps per track. Writing them now into
  //                       file...");
  //
  // // Write the sim hit track steps info to the event store
  // context.eventStore.add(m_cfg.geantTrackStepCollection,
  // std::move(trackSteps));
  return FW::ProcessCode::SUCCESS;
}
