// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

namespace ActsFatras {

/// Either return the global run manager or create a new dummy one.
///
/// @return Pointer to the run manager
G4RunManager* ensureGeant4RunManager();

/// Construct a dummy Geant4 detector.
class DummyDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  /// Destructor
  ~DummyDetectorConstruction() override;

  /// Interface method for Geant4.
  G4VPhysicalVolume* Construct() override;

 private:
  /// This method builds a dummy Geant4 detector.
  void dummyDetector();

  // Logical volume
  G4LogicalVolume* m_worldLog = nullptr;

  // Physical volume
  G4VPhysicalVolume* m_worldPhys = nullptr;
};
}  // namespace ActsFatras
