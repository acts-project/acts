// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <stdexcept>

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>
#include <globals.hh>

ActsExamples::SensitiveSurfaceMapper::SensitiveSurfaceMapper(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("No Acts::TrackingGeometry provided.");
  }
}

void ActsExamples::SensitiveSurfaceMapper::remapSensitiveNames(
    G4VPhysicalVolume* g4PhysicalVolume, Acts::Transform3 motherTransform,
    int& sCounter) const {
  auto g4LogicalVolume = g4PhysicalVolume->GetLogicalVolume();
  auto g4SensitiveDetector = g4LogicalVolume->GetSensitiveDetector();

  G4int nDaughters = g4LogicalVolume->GetNoDaughters();

  constexpr double convertLength = CLHEP::mm / Acts::UnitConstants::mm;

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

  if (nDaughters > 0) {
    // Step down to all daughters
    for (G4int id = 0; id < nDaughters; ++id) {
      remapSensitiveNames(g4LogicalVolume->GetDaughter(id), transform,
                          sCounter);
    }
    return;
  }

  std::string volumeName = g4LogicalVolume->GetName();
  std::string volumeMaterialName = g4LogicalVolume->GetMaterial()->GetName();
  if (g4SensitiveDetector != nullptr or
      std::find(m_cfg.materialMappings.begin(), m_cfg.materialMappings.end(),
                volumeMaterialName) != m_cfg.materialMappings.end() or
      std::find(m_cfg.volumeMappings.begin(), m_cfg.volumeMappings.end(),
                volumeName) != m_cfg.volumeMappings.end()) {
    // Find the associated ACTS object
    auto actsLayer = m_cfg.trackingGeometry->associatedLayer(
        Acts::GeometryContext(), g4AbsPosition);

    // Prepare the mapped surface
    const Acts::Surface* mappedSurface = nullptr;

    if (actsLayer != nullptr and actsLayer->surfaceArray() != nullptr) {
      auto actsSurfaces = actsLayer->surfaceArray()->at(g4AbsPosition);
      if (not actsSurfaces.empty()) {
        // Fast matching: search
        for (const auto& as : actsSurfaces) {
          if (as->center(Acts::GeometryContext()).isApprox(g4AbsPosition)) {
            mappedSurface = as;
            break;
          }
        }
      }
      if (mappedSurface == nullptr) {
        // Slow matching: Fallback, loop over all layer surfaces
        for (const auto& as : actsLayer->surfaceArray()->surfaces()) {
          if (as->center(Acts::GeometryContext()).isApprox(g4AbsPosition)) {
            mappedSurface = as;
            break;
          }
        }
      }
    }
    // A mapped surface was found, a new name will be set that
    // contains the GeometryID/
    if (mappedSurface != nullptr) {
      ++sCounter;
      std::string mappedVolumeName = SensitiveSurfaceMapper::mappingPrefix;
      mappedVolumeName += std::to_string(mappedSurface->geometryId().value());
      ACTS_VERBOSE("Found matching surface " << mappedSurface->geometryId()
                                             << " at position "
                                             << g4RelPosition.transpose());
      ACTS_VERBOSE("Remap: " << g4PhysicalVolume->GetName() << " -> "
                             << mappedVolumeName);
      g4PhysicalVolume->SetName(mappedVolumeName.c_str());
    } else {
      ACTS_VERBOSE("No mapping found for '"
                   << volumeName << "' with material '" << volumeMaterialName
                   << "' at position " << g4RelPosition.transpose());
    }
  } else {
    ACTS_VERBOSE(
        "Did not try mapping '"
        << g4PhysicalVolume->GetName() << "' at " << g4RelPosition.transpose()
        << " because g4SensitiveDetector is nullptr"
        << " and none of the provided volumes or materials were found");
  }
}
