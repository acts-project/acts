// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Plugins/Geant4/G4DetectorConstruction.hpp"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"

ActsFatras::G4DetectorConstruction::~G4DetectorConstruction() {
  delete (m_worldLog);
  delete (m_worldPhys);
}

void ActsFatras::G4DetectorConstruction::dummyDetector() {
  G4ThreeVector materialPosition(0., 0., 0.);

  // Create the world setup
  G4Box* worldBox = new G4Box("WorldBox", 25000., 25000., 25000.);

  // G4 material: vacuum setup
  G4Material* g4vacuum = G4Material::GetMaterial("Vacuum", false);
  if (!g4vacuum)
    g4vacuum =
        new G4Material("FatrasDummyVacuum", 1., 1.01 * CLHEP::g / CLHEP::mole,
                       CLHEP::universe_mean_density, kStateGas,
                       0.1 * CLHEP::kelvin, 1.e-19 * CLHEP::pascal);

  // Build the logical and physical volume
  m_worldLog =
      m_worldLog
          ? new (m_worldLog)
                G4LogicalVolume(worldBox, g4vacuum, "WorldLogical", 0, 0, 0)
          : new G4LogicalVolume(worldBox, g4vacuum, "WorldLogical", 0, 0, 0);
  m_worldPhys =
      m_worldPhys ? new (m_worldPhys)
                        G4PVPlacement(0, G4ThreeVector(0., 0., 0),
                                      "WorldPhysical", m_worldLog, 0, false, 0)
                  : new G4PVPlacement(0, materialPosition, "WorldPhysical",
                                      m_worldLog, 0, false, 0);
}

G4VPhysicalVolume* ActsFatras::G4DetectorConstruction::Construct() {
  // Construct the detector and return it
  dummyDetector();
  return m_worldPhys;
}
