// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <fstream>
#include <map>
#include <mutex>

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Plugins/Json/JsonGeometryConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"

// Convenience shorthand

namespace Acts {

/// @brief Material decorator from Json format
///
/// This reads in material maps for surfaces and volumes
/// from a json file
class JsonMaterialDecorator : public IMaterialDecorator {
 public:
  using SurfaceMaterialMap =
      std::map<GeometryID, std::shared_ptr<const ISurfaceMaterial>>;

  using VolumeMaterialMap =
      std::map<GeometryID, std::shared_ptr<const IVolumeMaterial>>;

  JsonMaterialDecorator(const JsonGeometryConverter::Config& rConfig,
                        const std::string& jFileName,
                        bool clearSurfaceMaterial = true,
                        bool clearVolumeMaterial = true)
      : m_readerConfig(rConfig),
        m_clearSurfaceMaterial(clearSurfaceMaterial),
        m_clearVolumeMaterial(clearVolumeMaterial) {
    // the material reader
    Acts::JsonGeometryConverter::Config jmConverterCfg("JsonGeometryConverter",
                                                       Logging::VERBOSE);
    Acts::JsonGeometryConverter jmConverter(jmConverterCfg);

    std::ifstream ifj(jFileName.c_str());
    nlohmann::json jin;
    ifj >> jin;

    auto maps = jmConverter.jsonToMaterialMaps(jin);
    m_surfaceMaterialMap = maps.first;
    m_volumeMaterialMap = maps.second;
  }

  /// Decorate a surface
  ///
  /// @param surface the non-cost surface that is decorated
  void decorate(Surface& surface) const final {
    // Clear the material if registered to do so
    if (m_clearSurfaceMaterial) {
      surface.assignSurfaceMaterial(nullptr);
    }
    // Try to find the surface in the map
    auto sMaterial = m_surfaceMaterialMap.find(surface.geoID());
    if (sMaterial != m_surfaceMaterialMap.end()) {
      surface.assignSurfaceMaterial(sMaterial->second);
    }
  }

  /// Decorate a TrackingVolume
  ///
  /// @param volume the non-cost volume that is decorated
  void decorate(TrackingVolume& volume) const final {
    // Clear the material if registered to do so
    if (m_clearVolumeMaterial) {
      volume.assignVolumeMaterial(nullptr);
    }
    // Try to find the volume in the map
    auto vMaterial = m_volumeMaterialMap.find(volume.geoID());
    if (vMaterial != m_volumeMaterialMap.end()) {
      volume.assignVolumeMaterial(vMaterial->second);
    }
  }

 private:
  JsonGeometryConverter::Config m_readerConfig;
  SurfaceMaterialMap m_surfaceMaterialMap;
  VolumeMaterialMap m_volumeMaterialMap;

  bool m_clearSurfaceMaterial{true};
  bool m_clearVolumeMaterial{true};
};
}  // namespace Acts
