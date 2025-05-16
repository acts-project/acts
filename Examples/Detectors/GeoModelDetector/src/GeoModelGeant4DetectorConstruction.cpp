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
#include <GeoModelKernel/throwExcept.h>
namespace ActsExamples {

GeoModelGeant4DetectorConstruction::GeoModelGeant4DetectorConstruction(const Config& cfg)
    : G4VUserDetectorConstruction(),
      m_options(cfg) {
  if (!m_options.geoModelTree.worldVolume) {
    throw std::invalid_argument(
        "GeoModelGeant4DetectorConstruction: "
        "GeoModel world volume is nullptr");
  }
}
void GeoModelGeant4DetectorConstruction::markSensitiveVols(G4VPhysicalVolume* g4Vol,
                                                           std::unordered_set<G4LogicalVolume*>& processed) const {
  //// Check that the volume has not yet been processed
  G4LogicalVolume* logVol = g4Vol->GetLogicalVolume();
  if (!processed.insert(logVol).second) {
      return;
  }
  const auto& volName = logVol->GetName();
  if (std::ranges::any_of(m_options.sensitiveVols, [&volName](const std::string& sensitiveToken){
      return volName.find(sensitiveToken) != std::string::npos; })) {
      /// The volume is sensitive -> don't traverse the daughters further  
      g4Vol->SetName(std::string{mappingPrefix} + g4Vol->GetName());
      return;
  }
  for (std::size_t daught=0 ; daught < logVol->GetNoDaughters(); ++daught) {
      markSensitiveVols(logVol->GetDaughter(daught), processed);
  }
}
G4VPhysicalVolume* GeoModelGeant4DetectorConstruction::Construct() {
  if (m_g4World == nullptr) {
    ExtParameterisedVolumeBuilder builder(m_options.geoModelTree.worldVolumeName);
    G4LogicalVolume* g4WorldLog = builder.Build(m_options.geoModelTree.worldVolume);
    
    std::unordered_set<G4LogicalVolume*> processedVols{};

    m_g4World =
        new G4PVPlacement(nullptr, G4ThreeVector(), g4WorldLog,
                          m_options.geoModelTree.worldVolumeName, nullptr, false, 0);

    markSensitiveVols(m_g4World, processedVols);
    // Create regions
    for (const auto& regionCreator : m_options.regionCreators) {
      regionCreator->buildRegion();
    }
  }
  return m_g4World;
}

}  // namespace ActsExamples
