// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <G4VUserDetectorConstruction.hh>

namespace dd4hep {
class Detector;
}
namespace ActsExamples {

/// Construct the Geant4 detector from a DD4hep description.
class DDG4DetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  DDG4DetectorConstruction(dd4hep::Detector& detector);

  /// Convert the stored DD4hep detector to a Geant4 description.
  ///
  /// Transfers ownership of the created object as all volumes (including world)
  /// are deleted in ~G4PhysicalVolumeStore().
  ///
  /// @note for facilitating configuration within the ACTS framework the world
  /// volume is cached
  G4VPhysicalVolume* Construct() final;

 private:
  /// The DD4hep detector instrance
  dd4hep::Detector& m_detector;
  /// The world volume
  G4VPhysicalVolume* m_world = nullptr;
};

}  // namespace ActsExamples
