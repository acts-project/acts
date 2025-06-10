// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"

#include <utility>

#include <G4ParticleTypes.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>
#include <G4UnitsTable.hh>

namespace ActsExamples::Geant4 {

MaterialPhysicsList::MaterialPhysicsList(
    std::unique_ptr<const Acts::Logger> logger)
    : G4VUserPhysicsList(), m_logger(std::move(logger)) {
  defaultCutValue = 1.0 * CLHEP::cm;
}

void MaterialPhysicsList::ConstructParticle() {
  ACTS_DEBUG("Construct Geantinos and Charged Geantinos.");
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
}

void MaterialPhysicsList::ConstructProcess() {
  ACTS_DEBUG("Adding Transport as single supperted Process.");
  AddTransportation();
}

void MaterialPhysicsList::SetCuts() {
  SetCutsWithDefault();

  if (verboseLevel > 0) {
    DumpCutValuesTable();
  }
}

}  // namespace ActsExamples::Geant4
