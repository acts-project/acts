// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TelescopeDetector/TelescopeG4DetectorConstruction.hpp"

#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/Geant4/RegionCreator.hpp"
#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"

#include <algorithm>
#include <cmath>

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

namespace ActsExamples {

TelescopeG4DetectorConstruction::TelescopeG4DetectorConstruction(
    const TelescopeDetector::Config& cfg,
    const Geant4ConstructionOptions& options)
    : m_cfg(cfg), m_options(options) {
  throw_assert(cfg.surfaceType == static_cast<int>(TelescopeSurfaceType::Plane),
               "only plane is supported right now");
}

G4VPhysicalVolume* TelescopeG4DetectorConstruction::Construct() {
  if (m_world != nullptr) {
    return m_world;
  }

  // Option to switch on/off checking of volumes overlaps
  constexpr G4bool checkOverlaps = true;

  // Keeps nested volumes from sharing surfaces
  constexpr G4double margin = 1 * mm;

  // The positions need not be sorted
  const auto [minPosition, maxPosition] = std::ranges::minmax(m_cfg.positions);

  // Extent of the layer stack along the telescope axis
  const G4double stackCenter = (minPosition + maxPosition) * 0.5 * mm;
  const G4double stackLength = (maxPosition - minPosition) * mm;

  // `bounds` are half lengths, matching the `Acts::RectangleBounds` that
  // `buildTelescopeDetector` gives the sensitive surfaces
  const G4double layerHalfX = m_cfg.bounds[0] * mm;
  const G4double layerHalfY = m_cfg.bounds[1] * mm;
  const G4double layerHalfZ = m_cfg.thickness * 0.5 * mm;

  const G4double envHalfX = layerHalfX + margin;
  const G4double envHalfY = layerHalfY + margin;
  const G4double envHalfZ = stackLength * 0.5 + layerHalfZ + margin;

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // `SensitiveSurfaceMapper` keys on the silicon material name
  G4Material* galactic = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* silicon = G4Material::GetMaterial("Silicon", false);
  if (silicon == nullptr) {
    silicon =
        new G4Material("Silicon", 14, 28.0855 * g / mole, 2.329 * g / cm3);
  }

  // Orientation of the envelope in the world, assuming binValue is AxisX,
  // AxisY or AxisZ. `G4PVPlacement` takes the rotation of the mother relative
  // to the daughter frame, i.e. the inverse of that orientation.
  G4RotationMatrix* rotation = nullptr;
  if (static_cast<Acts::AxisDirection>(m_cfg.binValue) ==
      Acts::AxisDirection::AxisX) {
    rotation = new G4RotationMatrix({0, 0, 1}, {0, 1, 0}, {-1, 0, 0});
  } else if (static_cast<Acts::AxisDirection>(m_cfg.binValue) ==
             Acts::AxisDirection::AxisY) {
    rotation = new G4RotationMatrix({1, 0, 0}, {0, 0, 1}, {0, -1, 0});
  }

  // The envelope center in the world frame
  G4ThreeVector envCenter(m_cfg.offsets[0] * mm, m_cfg.offsets[1] * mm,
                          stackCenter);
  if (rotation != nullptr) {
    envCenter = rotation->inverse() * envCenter;
  }

  // The world has to contain the envelope in any orientation
  const G4double worldHalfSize =
      envCenter.mag() +
      std::sqrt(envHalfX * envHalfX + envHalfY * envHalfY +
                envHalfZ * envHalfZ) +
      margin;

  // World
  //
  G4Box* solidWorld =
      new G4Box("World Solid", worldHalfSize, worldHalfSize, worldHalfSize);

  G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld, galactic, "World Logic");

  m_world = new G4PVPlacement(nullptr,          // no rotation
                              G4ThreeVector(),  // position
                              logicWorld,       // its logical volume
                              "World Phys",     // its name
                              nullptr,          // its mother volume
                              false,            // no boolean operation
                              0,                // copy number
                              checkOverlaps);   // overlaps checking

  // Envelope
  //
  G4Box* solidEnv = new G4Box("Envelope Solid",               // its name
                              envHalfX, envHalfY, envHalfZ);  // its size

  G4LogicalVolume* logicEnv =
      new G4LogicalVolume(solidEnv,           // its solid
                          galactic,           // its material
                          "Envelope Logic");  // its name

  G4VPhysicalVolume* physEnv =
      new G4PVPlacement(rotation,         // rotation
                        envCenter,        // at the center of the stack
                        logicEnv,         // its logical volume
                        "Envelope Phys",  // its name
                        logicWorld,       // its mother volume
                        false,            // no boolean operation
                        0,                // copy number
                        checkOverlaps);   // overlaps checking

  // Layer
  //

  G4Box* solidLayer =
      new G4Box("Layer Solid", layerHalfX, layerHalfY, layerHalfZ);

  G4LogicalVolume* logicLayer = new G4LogicalVolume(solidLayer,  // its solid
                                                    silicon,     // its material
                                                    "Layer Logic");  // its name

  for (std::size_t i = 0; i < m_cfg.positions.size(); ++i) {
    new G4PVPlacement(
        nullptr,  // no rotation
        G4ThreeVector(0, 0, m_cfg.positions[i] * mm - stackCenter),  // position
        "Layer #" + std::to_string(i) + " Phys",                     // its name
        logicLayer,             // its logical volume
        physEnv,                // its mother volume
        false,                  // no boolean operation
        static_cast<G4int>(i),  // copy number
        checkOverlaps);         // overlaps checking
  }

  // Create regions
  for (const auto& regionCreator : m_options.regionCreators) {
    regionCreator->buildRegion();
  }

  return m_world;
}

}  // namespace ActsExamples
