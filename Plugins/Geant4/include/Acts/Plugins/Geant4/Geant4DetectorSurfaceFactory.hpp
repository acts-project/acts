// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <tuple>
#include <vector>

#include "G4Transform3D.hh"

class G4VPhysicalVolume;

namespace Acts {

class Geant4DetectorElement;

/// A factory to convert Geant4 physical volumes
/// into Geant4 detector elements
///
class Geant4DetectorSurfaceFactory {
 public:
  /// Nested configuration struct that holds
  /// global lifetime configuration
  struct Config {};

  // Collect the senstive surfaces
  using Geant4SensitiveSurface =
      std::tuple<std::shared_ptr<Geant4DetectorElement>,
                 std::shared_ptr<Surface>>;

  // Collect the passive surfaces
  using Geant4PassiveSurface = std::shared_ptr<Surface>;

  /// Nested cache that records the conversion status
  struct Cache {
    /// The created detector elements - for the detector store
    std::vector<Geant4SensitiveSurface> sensitiveSurfaces;
    /// The created non-const surfaces - for further processing,
    std::vector<Geant4PassiveSurface> passiveSurfaces;
    /// matching and conversion statistics: volumes
    std::size_t matchedG4Volumes = 0;
    /// matching and conversion statistics: surfaces
    std::size_t convertedSurfaces = 0;
    /// matching and conversion statistics: materials
    std::size_t convertedMaterials = 0;
  };

  /// Nested option struct that allows per call changable configuration
  struct Options {
    /// Convert the length scale
    ActsScalar scaleConversion = 1.;
    /// Convert the material
    bool convertMaterial = true;
    /// Converted material thickness (< 0 indicates keeping original thickness)
    ActsScalar convertedMaterialThickness = -1;
    /// A selector for sensitive surfaces
    std::shared_ptr<IGeant4PhysicalVolumeSelector> sensitiveSurfaceSelector =
        nullptr;
    /// A selector for passive surfaces
    std::shared_ptr<IGeant4PhysicalVolumeSelector> passiveSurfaceSelector =
        nullptr;
  };

  /// The Geant4 detector element factory
  ///
  /// @param cfg the configuration struct
  Geant4DetectorSurfaceFactory() = default;

  /// Construction method of the detector elements
  ///
  /// @param cache [in,out] into which the Elements are filled
  /// @param g4ToGlobal the transformation to global
  /// @param g4PhyVol the current physical volume
  /// @param options the factory creation option
  ///
  void construct(Cache& cache, const G4Transform3D& g4ToGlobal,
                 const G4VPhysicalVolume& g4PhysVol, const Options& option);
};
}  // namespace Acts
