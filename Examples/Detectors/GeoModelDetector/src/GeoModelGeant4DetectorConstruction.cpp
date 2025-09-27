// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelDetector/GeoModelGeant4DetectorConstruction.hpp"

#include "ActsExamples/Geant4/Geant4ConstructionOptions.hpp"
#include "ActsExamples/Geant4/RegionCreator.hpp"

#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4ThreeVector.hh>
#include <G4VPhysicalVolume.hh>
#include <GeoModel2G4/ExtParameterisedVolumeBuilder.h>
#include <GeoModelKernel/GeoFullPhysVol.h>

namespace ActsExamples {

GeoModelGeant4DetectorConstruction::GeoModelGeant4DetectorConstruction(
    const ActsPlugins::GeoModelTree& geoModelTree,
    const Geant4ConstructionOptions& options)
    : G4VUserDetectorConstruction(),
      m_geoModelTree(geoModelTree),
      m_options(options) {
  if (geoModelTree.worldVolume == nullptr) {
    throw std::invalid_argument(
        "GeoModelGeant4DetectorConstruction: "
        "GeoModel world volume is nullptr");
  }
}

G4VPhysicalVolume* GeoModelGeant4DetectorConstruction::Construct() {
  if (m_g4World == nullptr) {
    ExtParameterisedVolumeBuilder builder(m_geoModelTree.worldVolumeName);
    G4LogicalVolume* g4WorldLog = builder.Build(m_geoModelTree.worldVolume);
    m_g4World =
        new G4PVPlacement(nullptr, G4ThreeVector(), g4WorldLog,
                          m_geoModelTree.worldVolumeName, nullptr, false, 0);

    // Create regions
    for (const auto& regionCreator : m_options.regionCreators) {
      regionCreator->buildRegion();
    }
  }
  return m_g4World;
}

}  // namespace ActsExamples
