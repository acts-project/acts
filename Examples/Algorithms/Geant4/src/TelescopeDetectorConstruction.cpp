#include "ActsExamples/TelescopeDetector/TelescopeDetectorConstruction.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "globals.hh"

ActsExamples::Telescope::TelescopeG4DetectorConstruction::
    TelescopeG4DetectorConstruction(TelescopeDetector::Config cfg)
    : m_cfg(cfg) {
  throw_assert(
      cfg.surfaceType == static_cast<int>(Telescope::TelescopeSurfaceType::Plane),
      "only plan is supported right now");
}

G4VPhysicalVolume*
ActsExamples::Telescope::TelescopeG4DetectorConstruction::Construct() {
  G4double length = m_cfg.positions.back() - m_cfg.positions.front();

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeX = m_cfg.bounds[0] + 5. * mm;
  G4double env_sizeY = m_cfg.bounds[1] + 5. * mm;
  G4double env_sizeZ = length + 10. * mm;

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // Materials
  G4Material* galactic = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");

  //
  // World
  //
  G4double worldSizeX = 2 * m_cfg.bounds[0];
  G4double worldSizeY = 2 * m_cfg.bounds[1];
  G4double worldSizeZ = length;

  // Construct the rotation
  // This assumes the binValue is binX, binY or binZ. No reset is necessary in
  // case of binZ
  G4RotationMatrix* rotation = nullptr;
  if (m_cfg.binValue == Acts::BinningValue::binX) {
    rotation = new G4RotationMatrix({0, 0, -1}, {0, 1, 0}, {1, 0, 0});
  } else if (m_cfg.binValue == Acts::BinningValue::binY) {
    rotation = new G4RotationMatrix({1, 0, 0}, {0, 0, -1}, {0, 1, 0});
  }

  G4Box* solidWorld =
      new G4Box("World", 0.5 * worldSizeX, 0.5 * worldSizeY, 0.5 * worldSizeZ);
  G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld, galactic, "World");

  G4VPhysicalVolume* physWorld =
      new G4PVPlacement(rotation,  // rotation
                        G4ThreeVector(m_cfg.offsets[0], m_cfg.offsets[1],
                                      length * 0.5),  // position
                        logicWorld,                   // its logical volume
                        "World",                      // its name
                        nullptr,                      // its mother  volume
                        false,                        // no boolean operation
                        0,                            // copy number
                        checkOverlaps);               // overlaps checking

  //
  // Envelope
  //
  G4Box* solidEnv =
      new G4Box("Envelope",                                          // its name
                0.5 * env_sizeX, 0.5 * env_sizeY, 0.5 * env_sizeZ);  // its size

  G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv,     // its solid
                                                  galactic,     // its material
                                                  "Envelope");  // its name

  new G4PVPlacement(nullptr,          // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicEnv,         // its logical volume
                    "Envelope",       // its name
                    logicWorld,       // its mother  volume
                    false,            // no boolean operation
                    0,                // copy number
                    checkOverlaps);   // overlaps checking

  for (std::size_t i = 0; i < m_cfg.positions.size(); ++i) {
    G4Box* surface = new G4Box("Surface #" + std::to_string(i), m_cfg.bounds[0],
                               m_cfg.bounds[1], m_cfg.thickness);

    G4LogicalVolume* locigShape =
        new G4LogicalVolume(surface,                         // its solid
                            silicon,                         // its material
                            "Shape #" + std::to_string(i));  // its name

    new G4PVPlacement(
        nullptr,                                       // no rotation
        G4ThreeVector(0, 0, m_cfg.positions[i] * mm),  // at position
        locigShape,                                    // its logical volume
        "Shape #" + std::to_string(i),                 // its name
        logicEnv,                                      // its mother  volume
        false,                                         // no boolean operation
        0,                                             // copy number
        checkOverlaps);                                // overlaps checking
  }

  return physWorld;
}
