// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

// see examples
// http://svn.code.sf.net/p/vgm/code/tags/v4-3/examples/E02/include/B4DetectorConstruction.hh

class G4VPhysicalVolume;
class TGeoNode;

namespace FW {
namespace Geant4 {

/// @class MMDetectorConstruction
///
/// @brief constructing the detector in Geant4
/// this can be used with GDML and TGeo iput
///
/// @TODO implement it with conversion from TGeo

class MMDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  /// Constructor
  MMDetectorConstruction();

  /// Construct method
  /// @return the world volume as G4VPhysicalVolume
  G4VPhysicalVolume* Construct() final override;

  /// Set the world TGeoNode to be transformed into
  /// a Geant4 geometry
  /// @param tgeoNode is the world not in ROOT::TGeo
  void setTGeoGeometry(TGeoNode* tgeoNode);

  /// Set the path and name to the GDML file
  /// @param gdmlFile is the path+name of the GDML file
  void setGdmlInput(std::string gdmlFile);

 private:
  TGeoNode* m_tgeoNode;
  std::string* m_gdmlFile;
};
}  // namespace Geant4
}  // namespace FW
