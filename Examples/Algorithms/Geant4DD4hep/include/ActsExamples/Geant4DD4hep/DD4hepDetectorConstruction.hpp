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

namespace dd4hep {
class Detector;
}
namespace ActsExamples {

/// Construct the Geant4 detector from a DD4hep description.
class DD4hepDetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  DD4hepDetectorConstruction(dd4hep::Detector& detector);

  /// Convert the stored DD4hep detector to a Geant4 description.
  ///
  /// Transfers ownership of the created object as all volumes (including world)
  /// are deleted in ~G4PhysicalVolumeStore().
  G4VPhysicalVolume* Construct() final override;

 private:
  dd4hep::Detector& m_detector;
};

class DD4hepDetectorConstructionFactory : public G4DetectorConstructionFactory {
 public:
  DD4hepDetectorConstructionFactory(dd4hep::Detector& detector);

  std::unique_ptr<G4VUserDetectorConstruction> operator()() const override;

 private:
  dd4hep::Detector& m_detector;
};

}  // namespace ActsExamples
