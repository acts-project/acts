// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/JsonGeometryConverter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Plugins/Json/MaterialJsonConverter.hpp"
#include "Acts/Plugins/Json/TrackingGeometryJsonConverter.hpp"

#include <map>

Acts::JsonGeometryConverter::JsonGeometryConverter(
    const Acts::JsonGeometryConverter::Config& cfg)
    : m_cfg(std::move(cfg)),
      m_volumeMaterialConverter(m_volumeName),
      m_volumeConverter(m_volumeName),
      m_surfaceMaterialConverter(m_surfaceName),
      m_surfaceConverter(m_surfaceName) {
  // Validate the configuration
  if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  }
}

/// Convert method
///
nlohmann::json Acts::JsonGeometryConverter::materialMapsToJson(
    const DetectorMaterialMaps& maps) {
  VolumeMaterialMap volumeMap = maps.second;
  std::vector<std::pair<GeometryIdentifier, const IVolumeMaterial*>>
      mapVolumeInit;
  for (auto it = volumeMap.begin(); it != volumeMap.end(); it++) {
    mapVolumeInit.push_back({it->first, it->second.get()});
  }
  GeometryHierarchyMap<const IVolumeMaterial*> HierarchyVolumeMap(
      mapVolumeInit);
  nlohmann::json materialVolume =
      m_volumeMaterialConverter.toJson(HierarchyVolumeMap);
  SurfaceMaterialMap surfaceMap = maps.first;
  std::vector<std::pair<GeometryIdentifier, const ISurfaceMaterial*>>
      mapSurfaceInit;
  for (auto it = surfaceMap.begin(); it != surfaceMap.end(); it++) {
    mapSurfaceInit.push_back({it->first, it->second.get()});
  }
  GeometryHierarchyMap<const ISurfaceMaterial*> HierarchySurfaceMap(
      mapSurfaceInit);
  nlohmann::json materialSurface =
      m_surfaceMaterialConverter.toJson(HierarchySurfaceMap);
  nlohmann::json materialMap;
  materialMap["1.Volumes"] = materialVolume;
  materialMap["2.Surfaces"] = materialSurface;
  return materialMap;
}

Acts::JsonGeometryConverter::DetectorMaterialMaps
Acts::JsonGeometryConverter::jsonToMaterialMaps(
    const nlohmann::json& materialmap) {
  nlohmann::json materialVolume = materialmap["1.Volumes"];
  GeometryHierarchyMap<const IVolumeMaterial*> HierarchyVolumeMap =
      m_volumeMaterialConverter.fromJson(materialVolume);
  VolumeMaterialMap volumeMap;
  for (size_t i = 0; i < HierarchyVolumeMap.size(); i++) {
    std::shared_ptr<const IVolumeMaterial> volumePointer(
        HierarchyVolumeMap.valueAt(i));
    volumeMap.insert({HierarchyVolumeMap.idAt(i), std::move(volumePointer)});
  }
  nlohmann::json materialSurface = materialmap["2.Surfaces"];
  GeometryHierarchyMap<const ISurfaceMaterial*> HierarchySurfaceMap =
      m_surfaceMaterialConverter.fromJson(materialSurface);
  SurfaceMaterialMap surfaceMap;
  for (size_t i = 0; i < HierarchySurfaceMap.size(); i++) {
    std::shared_ptr<const ISurfaceMaterial> surfacePointer(
        HierarchySurfaceMap.valueAt(i));
    surfaceMap.insert({HierarchySurfaceMap.idAt(i), std::move(surfacePointer)});
  }

  Acts::JsonGeometryConverter::DetectorMaterialMaps maps = {surfaceMap,
                                                            volumeMap};

  // Return the filled maps
  return maps;
}

nlohmann::json Acts::JsonGeometryConverter::trackingGeometryToJson(
    const Acts::TrackingGeometry& tGeometry) {
  std::vector<std::pair<GeometryIdentifier, const TrackingVolume*>>
      volumeHierarchy;
  std::vector<std::pair<GeometryIdentifier, const Surface*>> surfaceHierarchy;
  convertToHierarchy(volumeHierarchy, surfaceHierarchy,
                     tGeometry.highestTrackingVolume());
  GeometryHierarchyMap<const TrackingVolume*> HierarchyVolumeMap(
      volumeHierarchy);
  nlohmann::json jsonVolumes = m_volumeConverter.toJson(HierarchyVolumeMap);
  GeometryHierarchyMap<const Surface*> HierarchySurfaceMap(surfaceHierarchy);
  nlohmann::json jsonSurfaces = m_surfaceConverter.toJson(HierarchySurfaceMap);
  nlohmann::json hierarchyMap;
  hierarchyMap["1.Volumes"] = jsonVolumes;
  hierarchyMap["2.Surfaces"] = jsonSurfaces;
  return hierarchyMap;
}

void Acts::JsonGeometryConverter::convertToHierarchy(
    std::vector<std::pair<GeometryIdentifier, const TrackingVolume*>>&
        volumeHierarchy,
    std::vector<std::pair<GeometryIdentifier, const Surface*>>&
        surfaceHierarchy,
    const Acts::TrackingVolume* tVolume) {
  if ((tVolume->volumeMaterial() != nullptr ||
       m_cfg.processNonMaterial == true) &&
      m_cfg.processVolumes == true) {
    volumeHierarchy.push_back({tVolume->geometryId(), tVolume});
  }
  // there are confined volumes
  if (tVolume->confinedVolumes() != nullptr) {
    // get through the volumes
    auto& volumes = tVolume->confinedVolumes()->arrayObjects();
    // loop over the volumes
    for (auto& vol : volumes) {
      // recursive call
      convertToHierarchy(volumeHierarchy, surfaceHierarchy, vol.get());
    }
  }
  // there are dense volumes
  if (m_cfg.processDenseVolumes && !tVolume->denseVolumes().empty()) {
    // loop over the volumes
    for (auto& vol : tVolume->denseVolumes()) {
      // recursive call
      convertToHierarchy(volumeHierarchy, surfaceHierarchy, vol.get());
    }
  }
  if (tVolume->confinedLayers() != nullptr) {
    // get the layers
    auto& layers = tVolume->confinedLayers()->arrayObjects();
    // loop over the layers
    for (auto& lay : layers) {
      if (m_cfg.processRepresenting == true) {
        auto& layRep = lay->surfaceRepresentation();
        if ((layRep.surfaceMaterial() != nullptr ||
             m_cfg.processNonMaterial == true) &&
            layRep.geometryId() != GeometryIdentifier()) {
          surfaceHierarchy.push_back({layRep.geometryId(), &layRep});
        }
      }
      if (lay->approachDescriptor() != nullptr &&
          m_cfg.processApproaches == true) {
        for (auto& asf : lay->approachDescriptor()->containedSurfaces()) {
          if (asf->surfaceMaterial() != nullptr ||
              m_cfg.processNonMaterial == true) {
            surfaceHierarchy.push_back({asf->geometryId(), asf});
          }
        }
      }
      if (lay->surfaceArray() != nullptr && m_cfg.processSensitives == true) {
        for (auto& ssf : lay->surfaceArray()->surfaces()) {
          if (ssf->surfaceMaterial() != nullptr ||
              m_cfg.processNonMaterial == true) {
            surfaceHierarchy.push_back({ssf->geometryId(), ssf});
          }
        }
      }
    }
  }
  // Let's finally check the boundaries
  for (auto& bsurf : tVolume->boundarySurfaces()) {
    // the surface representation
    auto& bssfRep = bsurf->surfaceRepresentation();
    if (bssfRep.geometryId().volume() == tVolume->geometryId().volume() &&
        m_cfg.processBoundaries == true) {
      if (bssfRep.surfaceMaterial() != nullptr ||
          m_cfg.processNonMaterial == true) {
        surfaceHierarchy.push_back({bssfRep.geometryId(), &bssfRep});
      }
    }
  }
}
