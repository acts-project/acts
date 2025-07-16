// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Geant4/Geant4ConstructionOptions.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace ActsExamples {

class TelescopeG4DetectorConstruction final
    : public G4VUserDetectorConstruction {
 public:
  TelescopeG4DetectorConstruction(const TelescopeDetector::Config& cfg,
                                  const Geant4ConstructionOptions& options);

  G4VPhysicalVolume* Construct() final;

 private:
  /// The configuration of the telescope detector
  TelescopeDetector::Config m_cfg;
  /// The Geant4 construction options
  Geant4ConstructionOptions m_options;
  /// The world volume
  G4VPhysicalVolume* m_world{};
};

}  // namespace ActsExamples
