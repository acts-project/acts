// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Geant4/MMDetectorConstruction.hpp"

#include "G4GDMLParser.hh"
#include "TGeoManager.h"

FW::Geant4::MMDetectorConstruction::MMDetectorConstruction()
    : G4VUserDetectorConstruction(), m_tgeoNode(nullptr), m_gdmlFile(nullptr) {}

G4VPhysicalVolume* FW::Geant4::MMDetectorConstruction::Construct() {
  if (m_tgeoNode) {
    // Import geometry from Root to VGM
    /*     RootGM::Factory rtFactory;
         rtFactory.SetDebug(1);
         rtFactory.Import(m_tgeoNode);

         // Export VGM geometry to Geant4
         Geant4GM::Factory g4Factory;
         g4Factory.SetDebug(1);
         rtFactory.Export(&g4Factory);
         G4VPhysicalVolume* world = g4Factory.World();
         return world;*/
    return nullptr;
  } else if (m_gdmlFile) {
    G4GDMLParser parser;
    parser.Read(*m_gdmlFile);
    return parser.GetWorldVolume();
  } else
    return nullptr;  // and error Message
}

void FW::Geant4::MMDetectorConstruction::setTGeoGeometry(TGeoNode* tgeoNode) {
  m_tgeoNode = tgeoNode;
}

void FW::Geant4::MMDetectorConstruction::setGdmlInput(std::string gdmlFile) {
  m_gdmlFile = new std::string(gdmlFile);
}
