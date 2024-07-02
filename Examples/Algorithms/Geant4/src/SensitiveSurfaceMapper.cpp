// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include "Acts/Definitions/Units.hpp"

#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <utility>

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>

ActsExamples::SensitiveSurfaceMapper::SensitiveSurfaceMapper(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

void ActsExamples::SensitiveSurfaceMapper::remapSensitiveNames(
    State& state, const Acts::GeometryContext& gctx,
    G4VPhysicalVolume* g4PhysicalVolume,
    const Acts::Transform3& motherTransform) const {
  // Make sure the unit conversion is correct
  constexpr double convertLength = CLHEP::mm / Acts::UnitConstants::mm;

  auto g4LogicalVolume = g4PhysicalVolume->GetLogicalVolume();
  auto g4SensitiveDetector = g4LogicalVolume->GetSensitiveDetector();

  // Get the transform of the G4 object
  auto g4Translation = g4PhysicalVolume->GetTranslation();
  auto g4Rotation = g4PhysicalVolume->GetRotation();
  Acts::Vector3 g4RelPosition(g4Translation[0] * convertLength,
                              g4Translation[1] * convertLength,
                              g4Translation[2] * convertLength);
  Acts::Translation3 translation(g4RelPosition);
  Acts::Transform3 transform;
  if (g4Rotation == nullptr) {
    transform = motherTransform * translation;
  } else {
    Acts::RotationMatrix3 rotation;
    rotation << g4Rotation->xx(), g4Rotation->yx(), g4Rotation->zx(),
        g4Rotation->xy(), g4Rotation->yy(), g4Rotation->zy(), g4Rotation->xz(),
        g4Rotation->yz(), g4Rotation->zz();
    transform = motherTransform * (translation * rotation);
  }
  Acts::Vector3 g4AbsPosition = transform * Acts::Vector3::Zero();

  if (G4int nDaughters = g4LogicalVolume->GetNoDaughters(); nDaughters > 0) {
    // Step down to all daughters
    for (G4int id = 0; id < nDaughters; ++id) {
      remapSensitiveNames(state, gctx, g4LogicalVolume->GetDaughter(id),
                          transform);
    }
    return;
  }

  std::string volumeName = g4LogicalVolume->GetName();
  std::string volumeMaterialName = g4LogicalVolume->GetMaterial()->GetName();

  bool isSensitive = g4SensitiveDetector != nullptr;
  bool isMappedMaterial =
      std::find(m_cfg.materialMappings.begin(), m_cfg.materialMappings.end(),
                volumeMaterialName) != m_cfg.materialMappings.end();
  bool isMappedVolume =
      std::find(m_cfg.volumeMappings.begin(), m_cfg.volumeMappings.end(),
                volumeName) != m_cfg.volumeMappings.end();

  if (isSensitive || isMappedMaterial || isMappedVolume) {
    // Prepare the mapped surface
    const Acts::Surface* mappedSurface = nullptr;

    // Get the candidate surfaces
    auto candidateSurfaces = m_cfg.candidateSurfaces(gctx, g4AbsPosition);
    for (const auto& candidateSurface : candidateSurfaces) {
      // Check for center matching at the moment (needs to be extended)
      if (candidateSurface->center(gctx).isApprox(g4AbsPosition)) {
        mappedSurface = candidateSurface;
        break;
      }
    }

    // A mapped surface was found, a new name will be set that G4PhysVolume
    if (mappedSurface != nullptr) {
      ACTS_VERBOSE("Found matching surface " << mappedSurface->geometryId()
                                             << " at position "
                                             << g4RelPosition.transpose());
      // Check if the prefix is not yet assigned
      if (volumeName.find(mappingPrefix) == std::string::npos) {
        // Set the new name
        std::string mappedName = std::string(mappingPrefix) + volumeName;
        g4PhysicalVolume->SetName(mappedName);
      }
      // Insert into the multi-map
      state.g4VolumeToSurfaces.insert({g4PhysicalVolume, mappedSurface});
    } else {
      ACTS_VERBOSE("No mapping found for '"
                   << volumeName << "' with material '" << volumeMaterialName
                   << "' at position " << g4RelPosition.transpose());
    }
  } else {
    ACTS_VERBOSE("Did not try mapping '"
                 << g4PhysicalVolume->GetName() << "' at "
                 << g4RelPosition.transpose()
                 << " because g4SensitiveDetector (=" << g4SensitiveDetector
                 << ") is null and volume name (=" << volumeName
                 << ") and material name (=" << volumeMaterialName
                 << ") were not found");
  }
}
