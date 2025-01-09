// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Geant4/Geant4ConstructionOptions.hpp"

#include <memory>

#include <G4VUserDetectorConstruction.hh>

class G4VPhysicalVolume;

namespace dd4hep {
class Detector;
}  // namespace dd4hep

namespace ActsExamples {
class DD4hepDetector;

/// Construct the Geant4 detector from a DD4hep description.
class DDG4DetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  explicit DDG4DetectorConstruction(std::shared_ptr<dd4hep::Detector> detector,
                                    const Geant4ConstructionOptions& options);

  /// Convert the stored DD4hep detector to a Geant4 description.
  ///
  /// Transfers ownership of the created object as all volumes (including world)
  /// are deleted in ~G4PhysicalVolumeStore().
  ///
  /// @note for facilitating configuration within the ACTS framework the world
  /// volume is cached
  G4VPhysicalVolume* Construct() final;

 private:
  /// The DD4hep detector instance
  std::shared_ptr<dd4hep::Detector> m_detector;
  /// Construction options
  Geant4ConstructionOptions m_options;
  /// The world volume
  G4VPhysicalVolume* m_world = nullptr;
};

}  // namespace ActsExamples
