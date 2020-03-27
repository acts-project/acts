// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "DDG4/Geant4GeometryInfo.h"
#include "G4VUserDetectorConstruction.hh"

namespace dd4hep {
class Detector;
}
/// Temporary borrowed from FCCSW -> will be replaced later
/** @class GeoConstruction
 * DetectorDescription/DetDesServices/src/GeoConstruction.h GeoConstruction.h
 *
 *  Class to create Geant4 detector geometry from TGeo representation
 *  On demand (ie. when calling "Construct") the DD4hep geometry is converted
 *  to Geant4 with all volumes, assemblies, shapes, materials etc.
 *
 *  @author Markus Frank
 *  @author Anna Zaborowska
 */

namespace FW {

namespace DD4hepG4 {

class GeoConstruction : public G4VUserDetectorConstruction {
 public:
  /// Constructor
  GeoConstruction(dd4hep::Detector& lcdd);

  /// Geometry construction callback: Invoke the conversion to Geant4
  /// All volumes (including world) are deleted in ~G4PhysicalVolumeStore()
  G4VPhysicalVolume* Construct() final override;

 private:
  /// Reference to geometry object
  dd4hep::Detector& m_lcdd;
};
}  // namespace DD4hepG4
}  // namespace FW
