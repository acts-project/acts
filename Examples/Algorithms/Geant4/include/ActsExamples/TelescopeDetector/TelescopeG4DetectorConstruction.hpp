// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Geant4/DetectorConstructionFactory.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace ActsExamples {
namespace Telescope {

class TelescopeG4DetectorConstruction final
    : public G4VUserDetectorConstruction {
 public:
  TelescopeG4DetectorConstruction(const TelescopeDetector::Config& cfg);

  G4VPhysicalVolume* Construct() final;

 private:
  TelescopeDetector::Config m_cfg;

  G4VPhysicalVolume* m_world{};
};

class TelescopeG4DetectorConstructionFactory final
    : public DetectorConstructionFactory {
 public:
  TelescopeG4DetectorConstructionFactory(const TelescopeDetector::Config& cfg);

  std::unique_ptr<G4VUserDetectorConstruction> factorize() const override;

 private:
  TelescopeDetector::Config m_cfg;
};

}  // namespace Telescope
}  // namespace ActsExamples
