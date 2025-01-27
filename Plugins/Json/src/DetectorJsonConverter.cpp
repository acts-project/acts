// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/DetectorJsonConverter.hpp"

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Plugins/Json/DetectorVolumeFinderJsonConverter.hpp"
#include "Acts/Plugins/Json/DetectorVolumeJsonConverter.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"
#include "Acts/Plugins/Json/IndexedSurfacesJsonConverter.hpp"
#include "Acts/Plugins/Json/PortalJsonConverter.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <ctime>
#include <set>

nlohmann::json Acts::DetectorJsonConverter::toJson(
    const GeometryContext& gctx, const Experimental::Detector& detector,
    const Options& options) {
  // Get the time stamp
  std::time_t tt = 0;
  std::time(&tt);
  auto ti = std::localtime(&tt);

  nlohmann::json jDetector;

  std::size_t nSurfaces = 0;
  std::vector<const Experimental::Portal*> portals;

  for (const auto* volume : detector.volumes()) {
    nSurfaces += volume->surfaces().size();
    for (const auto& portal : volume->portals()) {
      if (std::find(portals.begin(), portals.end(), portal) == portals.end()) {
        portals.push_back(portal);
      }
    }
  }

  // Write the actual data objects
  nlohmann::json jData;
  jData["name"] = detector.name();

  // The portals are written first
  auto volumes = detector.volumes();
  nlohmann::json jPortals;

  for (const auto& portal : portals) {
    auto jPortal = PortalJsonConverter::toJson(
        gctx, *portal, volumes, options.volumeOptions.portalOptions);
    jPortals.push_back(jPortal);
  }
  jData["portals"] = jPortals;

  // The volumes are written next, with portals already defined, the volumes
  // will only index those
  nlohmann::json jVolumes;
  for (const auto& volume : volumes) {
    auto jVolume = DetectorVolumeJsonConverter::toJson(
        gctx, *volume, volumes, portals, options.volumeOptions);
    jVolumes.push_back(jVolume);
  }
  jData["volumes"] = jVolumes;

  jData["volume_finder"] = DetectorVolumeFinderJsonConverter::toJson(
      detector.detectorVolumeFinder());

  // Write the header
  nlohmann::json jHeader;
  jHeader["detector"] = detector.name();
  jHeader["type"] = "acts";
  jHeader["date"] = std::asctime(ti);
  jHeader["surface_count"] = nSurfaces;
  jHeader["portal_count"] = portals.size();
  jHeader["volume_count"] = detector.volumes().size();
  jDetector["header"] = jHeader;
  jDetector["data"] = jData;
  return jDetector;
}

nlohmann::json Acts::DetectorJsonConverter::toJsonDetray(
    const GeometryContext& gctx, const Experimental::Detector& detector,
    const Options& options) {
  // Get the time stamp
  std::time_t tt = 0;
  std::time(&tt);
  auto ti = std::localtime(&tt);

  nlohmann::json jFile;

  // (1) Detector section
  nlohmann::json jGeometry;
  nlohmann::json jGeometryData;
  nlohmann::json jGeometryHeader;
  std::size_t nSurfaces = 0;

  nlohmann::json jCommonHeader;
  jCommonHeader["detector"] = detector.name();
  jCommonHeader["date"] = std::asctime(ti);
  jCommonHeader["version"] = "detray - 0.44.0";
  jCommonHeader["tag"] = "geometry";

  auto volumes = detector.volumes();

  // Convert the volumes
  nlohmann::json jVolumes;
  for (const auto& volume : volumes) {
    auto jVolume = DetectorVolumeJsonConverter::toJsonDetray(
        gctx, *volume, volumes, options.volumeOptions);
    jVolumes.push_back(jVolume);
    if (jVolume.find("surfaces") != jVolume.end() &&
        jVolume["surfaces"].is_array()) {
      nSurfaces += jVolume["surfaces"].size();
    }
  }
  jGeometryData["volumes"] = jVolumes;
  jGeometryData["volume_grid"] = DetectorVolumeFinderJsonConverter::toJson(
      detector.detectorVolumeFinder(), true);

  // (1) Geometry
  // For detray, the number of surfaces and portals are collected
  // from the translated volumes
  jGeometryHeader["type"] = "detray";
  jGeometryHeader["common"] = jCommonHeader;
  jGeometryHeader["surface_count"] = nSurfaces;
  jGeometryHeader["volume_count"] = detector.volumes().size();
  jGeometry["header"] = jGeometryHeader;
  jGeometry["data"] = jGeometryData;
  jFile["geometry"] = jGeometry;

  // (2) surface grid section
  nlohmann::json jSurfaceGrids;
  nlohmann::json jSurfaceGridsData;
  nlohmann::json jSurfaceGridsCollection;
  nlohmann::json jSurfaceGridsHeader;
  for (const auto [iv, volume] : enumerate(volumes)) {
    // And its surface navigation delegates
    nlohmann::json jSurfacesDelegate = IndexedSurfacesJsonConverter::toJson(
        volume->surfaceCandidatesUpdater(), true);
    if (jSurfacesDelegate.is_null()) {
      continue;
    }
    // Patch axes for cylindrical grid surfaces, axes are swapped
    // at this point
    auto jAccLink = jSurfacesDelegate["acc_link"];
    std::size_t accLinkType = jAccLink["type"];
    if (accLinkType == 4u) {
      // Radial value to transfer phi to rphi
      std::vector<ActsScalar> bValues = volume->volumeBounds().values();
      ActsScalar rRef = 0.5 * (bValues[1] + bValues[0]);
      // Get the axes
      auto& jAxes = jSurfacesDelegate["axes"];
      // r*phi axis is the first one
      std::vector<ActsScalar> jAxesEdges = jAxes[0u]["edges"];
      std::for_each(jAxesEdges.begin(), jAxesEdges.end(),
                    [rRef](ActsScalar& phi) { phi *= rRef; });
      // Write back the patches axis edges
      jSurfacesDelegate["axes"][0u]["edges"] = jAxesEdges;
    }
    // Colplete the grid json for detray usage
    jSurfacesDelegate["volume_link"] = iv;
    // jSurfacesDelegate["acc_link"] =
    jSurfaceGridsCollection.push_back(jSurfacesDelegate);
  }
  jSurfaceGridsData["grids"] = jSurfaceGridsCollection;

  jCommonHeader["tag"] = "surface_grids";
  jSurfaceGridsHeader["common"] = jCommonHeader;
  jSurfaceGridsHeader["grid_count"] = jSurfaceGridsCollection.size();

  jSurfaceGrids["header"] = jSurfaceGridsHeader;
  jSurfaceGrids["data"] = jSurfaceGridsData;

  jFile["surface_grids"] = jSurfaceGrids;

  // (3) material section
  nlohmann::json jMaterial;

  jFile["material"] = jMaterial;

  return jFile;
}

std::shared_ptr<Acts::Experimental::Detector>
Acts::DetectorJsonConverter::fromJson(const GeometryContext& gctx,
                                      const nlohmann::json& jDetector) {
  // Read back all the data
  auto jData = jDetector["data"];
  auto jVolumes = jData["volumes"];
  auto jPortals = jData["portals"];
  const std::string name = jData["name"];

  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes;
  std::vector<std::shared_ptr<Experimental::Portal>> portals;

  for (const auto& jVolume : jVolumes) {
    auto volume = DetectorVolumeJsonConverter::fromJson(gctx, jVolume);
    volumes.push_back(volume);
  }

  for (const auto& jPortal : jPortals) {
    auto portal = PortalJsonConverter::fromJson(gctx, jPortal, volumes);
    portals.push_back(portal);
  }

  // Patch all the portals of the volumes
  for (auto [iv, v] : enumerate(volumes)) {
    // Get the Volumes and update the portals with the loaded ones
    auto jVolume = jVolumes[iv];
    std::vector<std::size_t> portalLinks = jVolume["portal_links"];
    for (auto [ip, ipl] : enumerate(portalLinks)) {
      auto portal = portals[ipl];
      v->updatePortal(portal, ip);
    }
  }
  return Experimental::Detector::makeShared(name, volumes,
                                            Experimental::tryRootVolumes());
}
