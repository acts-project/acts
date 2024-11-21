// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelDetector/GeoModelGeant4DetectorConstruction.hpp"

#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4ThreeVector.hh>
#include <G4VPhysicalVolume.hh>
#include <GeoModel2G4/ExtParameterisedVolumeBuilder.h>
#include <GeoModelKernel/GeoFullPhysVol.h>

namespace ActsExamples {

GeoModelGeant4DetectorConstruction::GeoModelGeant4DetectorConstruction(
    const Acts::GeoModelTree& geoModelTree,
    std::vector<std::shared_ptr<Geant4::RegionCreator>> regionCreators)
    : G4VUserDetectorConstruction(),
      m_geoModelTree(geoModelTree),
      m_regionCreators(std::move(regionCreators)) {
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
    for (const auto& regionCreator : m_regionCreators) {
      regionCreator->construct();
    }
  }
  return m_g4World;
}

GeoModelGeant4DetectorConstructionFactory::
    GeoModelGeant4DetectorConstructionFactory(
        const Acts::GeoModelTree& geoModelTree)
    : m_geoModelTree(geoModelTree) {}

std::unique_ptr<G4VUserDetectorConstruction>
GeoModelGeant4DetectorConstructionFactory::factorize(
    const std::vector<std::shared_ptr<Geant4::RegionCreator>>& regionCreators)
    const {
  return std::make_unique<GeoModelGeant4DetectorConstruction>(m_geoModelTree,
                                                              regionCreators);
}

}  // namespace ActsExamples
