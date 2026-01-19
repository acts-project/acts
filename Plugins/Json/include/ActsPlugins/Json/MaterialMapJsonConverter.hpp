// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/GeometryHierarchyMapJsonConverter.hpp"
#include "ActsPlugins/Json/ITrackingGeometryJsonDecorator.hpp"
#include "ActsPlugins/Json/IVolumeMaterialJsonDecorator.hpp"

#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{
class ISurfaceMaterial;
class ITrackingGeometryJsonDecorator;
class IVolumeMaterial;
class IVolumeMaterialJsonDecorator;
class Surface;
class TrackingGeometry;
class TrackingVolume;

/// A tuple containing a surface, its associated material, and geometry context.
/// @note Used for material mapping and JSON conversion.
using SurfaceAndMaterialWithContext =
    std::tuple<std::shared_ptr<const Acts::Surface>,
               std::shared_ptr<const Acts::ISurfaceMaterial>,
               Acts::GeometryContext>;
/// A pair containing a tracking volume and its associated material.
/// @note Used for material mapping and JSON conversion.
using TrackingVolumeAndMaterial =
    std::pair<const Acts::TrackingVolume*,
              std::shared_ptr<const Acts::IVolumeMaterial>>;

/// @class MaterialMapJsonConverter
///
/// @brief read the material from Json
class MaterialMapJsonConverter {
 public:
  /// @class Config
  /// Configuration of the Converter
  class Config {
   public:
    /// Default geometry context to extract surface transforms
    GeometryContext context = GeometryContext();

    /// Steering to handle sensitive data
    bool processSensitives = true;
    /// Steering to handle approach data
    bool processApproaches = true;
    /// Steering to handle representing data
    bool processRepresenting = true;
    /// Steering to handle boundary data
    bool processBoundaries = true;
    /// Steering to handle volume data
    bool processVolumes = true;
    /// Steering to handle volume data
    bool processDenseVolumes = false;
    /// Add proto material to all surfaces
    bool processNonMaterial = false;
  };

  /// Constructor
  ///
  /// @param config configuration struct for the reader
  /// @param level The log level
  MaterialMapJsonConverter(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~MaterialMapJsonConverter() = default;

  /// Convert a json material map to a TrackingGeometryMaterial
  ///
  /// @param materialmaps The json material
  /// @return Converted tracking geometry material from JSON
  TrackingGeometryMaterial jsonToMaterialMaps(
      const nlohmann::json& materialmaps);

  /// Convert a TrackingGeometryMaterial to json
  ///
  /// @param maps The material map collection
  /// @param decorator nullptr or a decorator to add extra attributes
  /// @return JSON representation of the material maps
  nlohmann::json materialMapsToJson(
      const TrackingGeometryMaterial& maps,
      const IVolumeMaterialJsonDecorator* decorator = nullptr);

  /// Convert a tracking geometry to json.
  /// Can be used to initialise the material mapping process.
  ///
  /// @param tGeometry is the tracking geometry
  /// @param decorator nullptr or a decorator to add extra attributes
  /// @return JSON representation of the tracking geometry
  nlohmann::json trackingGeometryToJson(
      const TrackingGeometry& tGeometry,
      const ITrackingGeometryJsonDecorator* decorator = nullptr);

  /// Go through a volume to find subvolume, layers and surfaces.
  /// Store volumes and surfaces in two vector used to initialised the geometry
  /// hierarchy.
  ///
  /// @param volumeHierarchy is a vector of volume to be filled
  /// @param surfaceHierarchy is a vector of surfaces to be filled
  /// @param tVolume is a volume
  void convertToHierarchy(
      std::vector<std::pair<GeometryIdentifier,
                            Acts::TrackingVolumeAndMaterial>>& volumeHierarchy,
      std::vector<
          std::pair<GeometryIdentifier, Acts::SurfaceAndMaterialWithContext>>&
          surfaceHierarchy,
      const Acts::TrackingVolume* tVolume);

 private:
  /// The config class
  Config m_cfg;

  /// The logger instance
  std::unique_ptr<const Logger> m_logger{nullptr};

  /// Name of the volume hierarchy
  std::string m_volumeName = "Material Volume Map";
  /// Geometry hierarchy writer for volume material.
  Acts::GeometryHierarchyMapJsonConverter<const IVolumeMaterial*,
                                          Acts::IVolumeMaterialJsonDecorator>
      m_volumeMaterialConverter;
  /// Geometry hierarchy writer for tracking volume.
  Acts::GeometryHierarchyMapJsonConverter<Acts::TrackingVolumeAndMaterial,
                                          Acts::ITrackingGeometryJsonDecorator>
      m_volumeConverter;

  /// Name of the surface hierarchy
  std::string m_surfaceName = "Material Surface Map";
  /// Geometry hierarchy writer for surface material.
  Acts::GeometryHierarchyMapJsonConverter<const ISurfaceMaterial*,
                                          Acts::IVolumeMaterialJsonDecorator>
      m_surfaceMaterialConverter;
  /// Geometry hierarchy writer for surface.
  Acts::GeometryHierarchyMapJsonConverter<Acts::SurfaceAndMaterialWithContext,
                                          Acts::ITrackingGeometryJsonDecorator>
      m_surfaceConverter;

  /// Private access to the logging instance
  const Logger& logger() const { return *m_logger; }
};

/// @}
}  // namespace Acts
