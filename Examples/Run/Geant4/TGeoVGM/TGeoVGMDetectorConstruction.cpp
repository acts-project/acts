// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TGeoVGMDetectorConstruction.hpp"

#include "Geant4GM/volumes/Factory.h"
#include "RootGM/volumes/Factory.h"
#include "TGeoManager.h"

using namespace ActsExamples;

TGeoVGMDetectorConstruction::TGeoVGMDetectorConstruction(std::string path)
    : G4VUserDetectorConstruction(), m_path(std::move(path)) {}

G4VPhysicalVolume* TGeoVGMDetectorConstruction::Construct() {
  // Import geometry from the root file
  new TGeoManager("TGeoManager", "Geant4 TGeo importer");
  gGeoManager->Import(m_path.c_str());

  // Import geometry from Root to VGM
  RootGM::Factory rtFactory;
  rtFactory.SetDebug(1);
  rtFactory.Import(gGeoManager->GetTopNode());

  // Export VGM geometry to Geant4
  //
  Geant4GM::Factory g4Factory;
  g4Factory.SetDebug(1);
  rtFactory.Export(&g4Factory);
  return g4Factory.World();
}
