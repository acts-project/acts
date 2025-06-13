// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cstddef>
#include <memory>
#include <tuple>
#include <vector>

namespace HepGeom {
class Transform3D;
}

class G4VPhysicalVolume;
using G4Transform3D = HepGeom::Transform3D;

namespace Acts {

class Geant4DetectorElement;
class IGeant4PhysicalVolumeSelector;
class Surface;

/// A factory to convert Geant4 physical volumes
/// into Geant4 detector elements
///
class Geant4DetectorSurfaceFactory {
 public:
  using ElementFactory = std::function<std::shared_ptr<Geant4DetectorElement>(
      std::shared_ptr<Surface>, const G4VPhysicalVolume&, const Transform3&,
      double)>;

  /// Nested configuration struct that holds
  /// global lifetime configuration
  struct Config {
    /// @cond
    /// The detector element factory with default implementation
    ElementFactory detectorElementFactory =
        [](std::shared_ptr<Surface> surface, const G4VPhysicalVolume& g4physVol,
           const Transform3& toGlobal, double thickness) {
          return std::make_shared<Geant4DetectorElement>(
              std::move(surface), g4physVol, toGlobal, thickness);
        };
    /// @endcond
  };

  // Collect the sensitive surfaces
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

  /// Nested option struct that allows per call changeable configuration
  struct Options {
    /// Convert the length scale
    double scaleConversion = 1.;
    /// Convert the material
    bool convertMaterial = false;
    /// Converted material thickness (< 0 indicates keeping original thickness)
    double convertedMaterialThickness = -1;
    /// A selector for sensitive surfaces
    std::shared_ptr<IGeant4PhysicalVolumeSelector> sensitiveSurfaceSelector =
        nullptr;
    /// A selector for passive surfaces
    std::shared_ptr<IGeant4PhysicalVolumeSelector> passiveSurfaceSelector =
        nullptr;
  };

  /// The Geant4 detector element factory
  Geant4DetectorSurfaceFactory() = delete;

  /// The Geant4 detector element factory with a configuration
  ///
  /// @param config the configuration of the factory
  /// @param mlogger a screen output logger
  explicit Geant4DetectorSurfaceFactory(
      const Config& config,
      std::unique_ptr<const Logger> mlogger =
          getDefaultLogger("Geant4DetectorSurfaceFactory", Acts::Logging::INFO))
      : m_config(config), m_logger(std::move(mlogger)) {}

  /// Construction method of the detector elements
  ///
  /// @param cache [in,out] into which the Elements are filled
  /// @param g4ToGlobal the transformation to global
  /// @param g4PhysVol the current physical volume
  /// @param option the factory creation option
  ///
  void construct(Cache& cache, const G4Transform3D& g4ToGlobal,
                 const G4VPhysicalVolume& g4PhysVol, const Options& option);

 private:
  /// The configuration of the factory
  Config m_config = Config{};

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
