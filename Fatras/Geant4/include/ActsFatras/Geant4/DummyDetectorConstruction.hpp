// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4RunManager.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

namespace ActsFatras {

/// @brief Convenience method assuring the existance of a G4RunManager
///
/// @return Pointer to the run manager
G4RunManager* makeDummyRunManager();

/// Construct a dummy Geant4 detector.
class DummyDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  /// Destructor
  ~DummyDetectorConstruction();

  /// @brief Interface method for Geant4
  G4VPhysicalVolume* Construct();

 private:
  /// @brief This method builds a dummy Geant4 detector
  void dummyDetector();

  // Logical volume
  G4LogicalVolume* m_worldLog = nullptr;

  // Physical volume
  G4VPhysicalVolume* m_worldPhys = nullptr;
};
}  // namespace ActsFatras