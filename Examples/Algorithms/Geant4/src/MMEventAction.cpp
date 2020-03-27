// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/MMEventAction.hpp"

#include <stdexcept>

#include "ACTFW/Plugins/Geant4/MMPrimaryGeneratorAction.hpp"
#include "ACTFW/Plugins/Geant4/MMSteppingAction.hpp"
#include "G4Event.hh"
#include "G4RunManager.hh"

FW::Geant4::MMEventAction* FW::Geant4::MMEventAction::fgInstance = nullptr;

FW::Geant4::MMEventAction* FW::Geant4::MMEventAction::Instance() {
  // Static acces function via G4RunManager
  return fgInstance;
}

FW::Geant4::MMEventAction::MMEventAction() : G4UserEventAction() {
  if (fgInstance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    fgInstance = this;
  }
}

FW::Geant4::MMEventAction::~MMEventAction() {
  fgInstance = nullptr;
}

void FW::Geant4::MMEventAction::BeginOfEventAction(const G4Event*) {
  // reset the collection of material steps
  MMSteppingAction::Instance()->Reset();
}

void FW::Geant4::MMEventAction::EndOfEventAction(const G4Event* event) {
  const auto* rawPos = event->GetPrimaryVertex();
  // access the initial direction of the track
  G4ThreeVector rawDir = MMPrimaryGeneratorAction::Instance()->direction();
  // create the RecordedMaterialTrack
  Acts::RecordedMaterialTrack mtrecord;
  mtrecord.first.first =
      Acts::Vector3D(rawPos->GetX0(), rawPos->GetY0(), rawPos->GetZ0());
  mtrecord.first.second = Acts::Vector3D(rawDir.x(), rawDir.y(), rawDir.z());
  mtrecord.second.materialInteractions =
      MMSteppingAction::Instance()->materialSteps();

  // write out the RecordedMaterialTrack of one event
  m_records.push_back(mtrecord);

  // write out the steps of one track in an event
  m_tracksteps.adopt_sequence(MMSteppingAction::Instance()->trackSteps());
}

void FW::Geant4::MMEventAction::Reset() {}
