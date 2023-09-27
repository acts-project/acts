// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"

#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

void Acts::Geant4DetectorSurfaceFactory::construct(
    Cache& cache, const G4Transform3D& g4ToGlobal,
    const G4VPhysicalVolume& g4PhysVol, const Options& option) {
  // Get Rotation and translation
  auto g4Translation = g4PhysVol.GetTranslation();
  auto g4Rotation = g4PhysVol.GetRotation();

  G4Transform3D g4Transform =
      (g4Rotation == nullptr)
          ? G4Transform3D(CLHEP::HepRotation(), g4Translation)
          : G4Transform3D(*g4Rotation, g4Translation);

  G4Transform3D newToGlobal = g4ToGlobal * g4Transform;

  // Get the logical volume
  auto g4LogicalVolume = g4PhysVol.GetLogicalVolume();
  std::size_t nDaughters = g4LogicalVolume->GetNoDaughters();
  for (std::size_t d = 0; d < nDaughters; ++d) {
    auto daughter = g4LogicalVolume->GetDaughter(d);
    construct(cache, newToGlobal, *daughter, option);
  }

  // Check if the volume is accepted by a sensitive or passive selector
  bool sensitive = option.sensitiveSurfaceSelector != nullptr and
                   option.sensitiveSurfaceSelector->select(g4PhysVol);
  bool passive = option.passiveSurfaceSelector != nullptr and
                 option.passiveSurfaceSelector->select(g4PhysVol);

  if (sensitive or passive) {
    // Conversion and selection code
    ++cache.matchedG4Volumes;

    // Attempt the conversion
    auto surface = Acts::Geant4PhysicalVolumeConverter{}.surface(
        g4PhysVol, Geant4AlgebraConverter{}.transform(newToGlobal),
        option.convertMaterial, option.convertedMaterialThickness);

    if (surface != nullptr) {
      ++cache.convertedSurfaces;
      // Count the material conversion
      if (surface->surfaceMaterial() != nullptr) {
        ++cache.convertedMaterials;
      }

      if (sensitive) {
        // empty gemetry context is fine as the transform was just passed down
        // without context before
        auto detectorElement = std::make_shared<Acts::Geant4DetectorElement>(
            surface, g4PhysVol, surface->transform({}), 0.1);
        surface->assignDetectorElement(*detectorElement);

        cache.sensitiveSurfaces.push_back(
            {std::move(detectorElement), std::move(surface)});
      } else {
        cache.passiveSurfaces.push_back(std::move(surface));
      }
    }
  }
}
