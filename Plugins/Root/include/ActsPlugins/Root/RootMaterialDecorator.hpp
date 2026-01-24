// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/RootMaterialMapIo.hpp"

#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <utility>

class TFile;

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// @class RootMaterialDecorator
///
/// @brief Read the collection of SurfaceMaterial & VolumeMaterial
class RootMaterialDecorator : public Acts::IMaterialDecorator {
 public:
  /// @class Config
  /// Configuration of the Reader
  class Config {
   public:
    /// Accessor config
    RootMaterialMapIo::Config accessorConfig;
    /// Accessor options
    RootMaterialMapIo::Options accessorOptions;
    /// The name of the output file
    std::string fileName = "material-maps.root";
  };

  /// Constructor
  ///
  /// @param config configuration struct for the reader
  /// @param level the logging level
  RootMaterialDecorator(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~RootMaterialDecorator() override;

  /// Decorate a surface
  ///
  /// @param surface the non-cost surface that is decorated
  void decorate(Acts::Surface& surface) const final {
    // Null out the material for this surface
    if (m_clearSurfaceMaterial) {
      surface.assignSurfaceMaterial(nullptr);
    }
    // Try to find the surface in the map
    auto sMaterial = m_surfaceMaterialMap.find(surface.geometryId());
    if (sMaterial != m_surfaceMaterialMap.end()) {
      surface.assignSurfaceMaterial(sMaterial->second);
    }
  }

  /// Decorate a TrackingVolume
  ///
  /// @param volume the non-cost volume that is decorated
  void decorate(Acts::TrackingVolume& volume) const final {
    // Null out the material for this volume
    if (m_clearSurfaceMaterial) {
      volume.assignVolumeMaterial(nullptr);
    }
    // Try to find the surface in the map
    auto vMaterial = m_volumeMaterialMap.find(volume.geometryId());
    if (vMaterial != m_volumeMaterialMap.end()) {
      volume.assignVolumeMaterial(vMaterial->second);
    }
  }

  /// Return the maps
  Acts::TrackingGeometryMaterial materialMaps() const {
    return {m_surfaceMaterialMap, m_volumeMaterialMap};
  }

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// The config class
  Config m_cfg;

  std::unique_ptr<const Acts::Logger> m_logger{nullptr};

  /// The input file
  TFile* m_inputFile{nullptr};

  /// Surface based material
  Acts::SurfaceMaterialMaps m_surfaceMaterialMap;

  /// Volume based material
  Acts::VolumeMaterialMaps m_volumeMaterialMap;

  bool m_clearSurfaceMaterial{true};

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }
};

/// @}
}  // namespace ActsPlugins
