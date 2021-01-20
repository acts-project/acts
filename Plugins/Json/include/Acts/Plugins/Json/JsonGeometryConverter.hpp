// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include "Acts/Plugins/Json/GeometryHierarchyMapJsonConverter.hpp"

#include <map>

#include <nlohmann/json.hpp>

namespace Acts {

/// @class JsonGeometryConverter
///
/// @brief read the material from Json
class JsonGeometryConverter {
 public:
  using SurfaceMaterialMap =
      std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;
  using VolumeMaterialMap =
      std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;
  using DetectorMaterialMaps = std::pair<SurfaceMaterialMap, VolumeMaterialMap>;

  /// @struct jsonKey
  ///
  /// @brief store in a single place the different key used for the material mapping 
  struct jsonKey {
    /// The name identification
    std::string namekey = "Name";
    /// The bin0 key
    std::string bin0key = "bin0";
    /// The bin1 key
    std::string bin1key = "bin1";
    /// The bin2 key
    std::string bin2key = "bin2";
    /// The local to global transformation key
    std::string transfokeys = "transformation";
    /// The type key -> proto, else
    std::string typekey = "type";
    /// The data key
    std::string datakey = "data";
    /// The geoid key
    std::string geometryidkey = "Geoid";
    /// The mapping key, add surface to mapping procedure if true
    std::string mapkey = "mapMaterial";
    /// The surface type key
    std::string surfacetypekey = "stype";
    /// The surface position key
    std::string surfacepositionkey = "sposition";
    /// The surface range key
    std::string surfacerangekey = "srange";
  };

  /// @class Config
  /// Configuration of the Converter
  class Config {
   public:
    /// The default logger
    std::shared_ptr<const Logger> logger;
    /// Default geometry context to extract surface tranforms
    GeometryContext context = GeometryContext();
    /// The name of the writer
    std::string name = "";

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

    /// Constructor
    ///
    /// @param lname Name of the writer tool
    /// @param lvl The output logging level
    Config(const std::string& lname = "JsonGeometryConverter",
           Logging::Level lvl = Logging::INFO)
        : logger(getDefaultLogger(lname, lvl)), name(lname) {}
  };

  /// Constructor
  ///
  /// @param cfg configuration struct for the reader
  JsonGeometryConverter(const Config& cfg);

  /// Destructor
  ~JsonGeometryConverter() = default;

  /// Convert a json material map to a DetectorMaterialMaps
  ///
  /// @param materialmaps The json material
  DetectorMaterialMaps jsonToMaterialMaps(const nlohmann::json& materialmaps);

  /// Convert a DetectorMaterialMaps to json
  ///
  /// @param maps The material map collection
  nlohmann::json materialMapsToJson(const DetectorMaterialMaps& maps);

  /// Convert a tracking geometry to json.
  /// Can be used to initialise the material mapping process.
  ///
  /// @param tGeometry is the tracking geometry
  nlohmann::json trackingGeometryToJson(const TrackingGeometry& tGeometry);

  /// Go through a volume to find subvolume, layers and surfaces.
  /// Store volumes and surfaces in two vector used to initialised the geometry hierachy.
  ///
  /// @param volumeHierarchy is a vector of volume to be filled
  /// @param surfaceHierarchy is a vector of surfaces to be filled
  /// @param tVolume is a volume 
  void convertToHierarchy(
      std::vector<std::pair<GeometryIdentifier, const TrackingVolume*>>&
          volumeHierarchy,
      std::vector<std::pair<GeometryIdentifier, const Surface*>>&
          surfaceHierarchy,
      const Acts::TrackingVolume* tVolume);

 private:

  /// The config class
  Config m_cfg;

  /// Name of the volume hierarchy
  std::string m_volumeName = "Material Volume Map";
  /// Geometry hierarchy writer for volume material.
  Acts::GeometryHierarchyMapJsonConverter<const IVolumeMaterial*> m_volumeMaterialConverter;
  /// Geometry hierarchy writer for tracking volume.
  Acts::GeometryHierarchyMapJsonConverter<const TrackingVolume*> m_volumeConverter;
  
  /// Name of the surface hierarchy
  std::string m_surfaceName = "Material Surface Map";
  /// Geometry hierarchy writer for surface material.
  Acts::GeometryHierarchyMapJsonConverter<const ISurfaceMaterial*> m_surfaceMaterialConverter;
/// Geometry hierarchy writer for surface.
  Acts::GeometryHierarchyMapJsonConverter<const Surface*> m_surfaceConverter;

  /// Private access to the logging instance
  const Logger& logger() const { return *m_cfg.logger; }
};

}  // namespace Acts
