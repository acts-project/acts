// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DetectorCommons/Geant4DetectorConstructionFactory.hpp"
#include "ActsExamples/Geant4/RegionCreator.hpp"

#include <memory>

#include <G4VUserDetectorConstruction.hh>

class G4VPhysicalVolume;

namespace ActsExamples {
class DD4hepDetector;

/// Construct the Geant4 detector from a DD4hep description.
class DDG4DetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  explicit DDG4DetectorConstruction(
      std::shared_ptr<DD4hepDetector> detector,
      std::vector<std::shared_ptr<Geant4::RegionCreator>> regionCreators);

  /// Convert the stored DD4hep detector to a Geant4 description.
  ///
  /// Transfers ownership of the created object as all volumes (including world)
  /// are deleted in ~G4PhysicalVolumeStore().
  ///
  /// @note for facilitating configuration within the ACTS framework the world
  /// volume is cached
  G4VPhysicalVolume* Construct() final;

 private:
  /// The Acts DD4hep detector instance
  std::shared_ptr<DD4hepDetector> m_detector;
  /// Region creators
  std::vector<std::shared_ptr<Geant4::RegionCreator>> m_regionCreators;
  /// The world volume
  G4VPhysicalVolume* m_world = nullptr;
};

class DDG4DetectorConstructionFactory final
    : public Geant4DetectorConstructionFactory {
 public:
  explicit DDG4DetectorConstructionFactory(
      std::shared_ptr<DD4hepDetector> detector);

  std::unique_ptr<G4VUserDetectorConstruction> factorize(
      const std::vector<std::shared_ptr<Geant4::RegionCreator>>& regionCreators)
      const override;

 private:
  std::shared_ptr<DD4hepDetector> m_detector;
};

}  // namespace ActsExamples
