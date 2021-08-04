// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Geant4/G4DetectorConstructionFactory.hpp"

#include <G4VUserDetectorConstruction.hh>
#include <globals.hh>

namespace ActsExamples {

/// Construct the Geant4 detector from a Gdml file.
class GdmlDetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  /// @param path is the path to the Gdml file
  GdmlDetectorConstruction(std::string path);

  /// Read the file and parse it to construct the Geant4 description.
  G4VPhysicalVolume* Construct() final override;

 private:
  std::string m_path;
};

class GdmlDetectorConstructionFactory : public G4DetectorConstructionFactory {
 public:
  GdmlDetectorConstructionFactory(const std::string& path);

  std::unique_ptr<G4VUserDetectorConstruction> operator()() const override;

 private:
  std::string m_path;
};

}  // namespace ActsExamples
