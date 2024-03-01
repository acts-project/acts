// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/DetectorVolumeJsonConverter.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdaters.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"
#include "Acts/Plugins/Json/IndexedSurfacesJsonConverter.hpp"
#include "Acts/Plugins/Json/PortalJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"
#include "Acts/Plugins/Json/VolumeBoundsJsonConverter.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <ctime>

namespace {

/// Find the position of the volume to point to
///
/// @param volume the volume to find
/// @param the collection of volumes
///
/// @note return -1 if not found, to be interpreted by the caller
int findVolume(
    const Acts::Experimental::DetectorVolume* volume,
    const std::vector<const Acts::Experimental::DetectorVolume*>& volumes) {
  auto candidate = std::find(volumes.begin(), volumes.end(), volume);
  if (candidate != volumes.end()) {
    return std::distance(volumes.begin(), candidate);
  }
  return -1;
}
}  // namespace

nlohmann::json Acts::DetectorVolumeJsonConverter::toJson(
    const GeometryContext& gctx, const Experimental::DetectorVolume& volume,
    const std::vector<const Experimental::DetectorVolume*>& detectorVolumes,
    const std::vector<const Experimental::Portal*>& portals,
    const Options& options) {
  nlohmann::json jVolume;
  jVolume["name"] = volume.name();
  jVolume["geometryId"] = volume.geometryId().volume();
  jVolume["transform"] = Transform3JsonConverter::toJson(
      volume.transform(gctx), options.transformOptions);
  jVolume["bounds"] = VolumeBoundsJsonConverter::toJson(volume.volumeBounds());
  // Write the surfaces
  nlohmann::json jSurfaces;
  std::for_each(
      volume.surfaces().begin(), volume.surfaces().end(), [&](const auto& s) {
        jSurfaces.push_back(
            SurfaceJsonConverter::toJson(gctx, *s, options.surfaceOptions));
      });
  jVolume["surfaces"] = jSurfaces;
  // And its surface navigation delegates
  nlohmann::json jSurfacesDelegate =
      IndexedSurfacesJsonConverter::toJson(volume.surfaceCandidatesUpdater());
  jVolume["surface_navigation"] = jSurfacesDelegate;

  // Write the sub volumes
  nlohmann::json jVolumes;
  std::for_each(
      volume.volumes().begin(), volume.volumes().end(), [&](const auto& v) {
        jVolumes.push_back(toJson(gctx, *v, detectorVolumes, portals, options));
      });
  jVolume["volumes"] = jVolumes;

  // Write the portals if pre-converted as link
  nlohmann::json jPortals;
  if (!portals.empty()) {
    for (const auto* p : volume.portals()) {
      auto it = std::find(portals.begin(), portals.end(), p);
      if (it != portals.end()) {
        jPortals.push_back(std::distance(portals.begin(), it));
      } else {
        throw std::runtime_error("Portal not found in the list of portals");
      }
    }
    jVolume["portal_links"] = jPortals;
  } else {
    for (const auto& p : volume.portals()) {
      nlohmann::json jPortal = PortalJsonConverter::toJson(
          gctx, *p, detectorVolumes, options.portalOptions);
      jPortals.push_back(jPortal);
    }
    jVolume["portals"] = jPortals;
  }
  return jVolume;
}

nlohmann::json Acts::DetectorVolumeJsonConverter::toJsonDetray(
    const GeometryContext& gctx, const Experimental::DetectorVolume& volume,
    const std::vector<const Experimental::DetectorVolume*>& detectorVolumes,
    const Options& options) {
  nlohmann::json jVolume;
  jVolume["name"] = volume.name();

  // Write the transform - path them with defaults
  jVolume["transform"] = Transform3JsonConverter::toJson(
      volume.transform(gctx), options.transformOptions);
  jVolume["bounds"] = VolumeBoundsJsonConverter::toJson(volume.volumeBounds());
  auto volumeBoundsType = volume.volumeBounds().type();
  if (volumeBoundsType == VolumeBounds::BoundsType::eCylinder) {
    jVolume["type"] = 0u;
  } else if (volumeBoundsType == VolumeBounds::BoundsType::eCuboid) {
    jVolume["type"] = 4u;
  } else {
    throw std::runtime_error("Unsupported volume bounds type");
  }

  // Get the index
  int vIndex = findVolume(&volume, detectorVolumes);
  jVolume["index"] = vIndex;

  std::size_t sIndex = 0;
  // Write the surfaces - patch bounds & augment with self links
  nlohmann::json jSurfaces;
  for (const auto& s : volume.surfaces()) {
    auto jSurface =
        SurfaceJsonConverter::toJsonDetray(gctx, *s, options.surfaceOptions);
    DetrayJsonHelper::addVolumeLink(jSurface["mask"], vIndex);
    jSurface["index_in_coll"] = sIndex++;
    jSurfaces.push_back(jSurface);
  }

  // Create the oriented surfaces, they could potentially be one-to-one
  // translated
  auto orientedSurfaces =
      volume.volumeBounds().orientedSurfaces(volume.transform(gctx));

  // Write the portals - they will end up in the surface container
  for (const auto& [ip, p] : enumerate(volume.portals())) {
    auto jPortalSurfaces =
        (toJsonDetray(gctx, *p, ip, volume, orientedSurfaces, detectorVolumes,
                      options.portalOptions));
    std::for_each(jPortalSurfaces.begin(), jPortalSurfaces.end(),
                  [&](auto& jSurface) {
                    jSurface["index_in_coll"] = sIndex++;
                    jSurfaces.push_back(jSurface);
                  });
  }
  jVolume["surfaces"] = jSurfaces;

  return jVolume;
}

std::shared_ptr<Acts::Experimental::DetectorVolume>
Acts::DetectorVolumeJsonConverter::fromJson(const GeometryContext& gctx,
                                            const nlohmann::json& jVolume) {
  std::string name = jVolume["name"];
  GeometryIdentifier geoId;
  geoId.setVolume(jVolume["geometryId"]);
  Transform3 transform =
      Transform3JsonConverter::fromJson(jVolume["transform"]);
  auto bounds = VolumeBoundsJsonConverter::fromJson(jVolume["bounds"]);

  auto jSurfaces = jVolume["surfaces"];
  auto jVolumes = jVolume["volumes"];

  // Some tooling
  auto portalGenerator = Experimental::defaultPortalGenerator();

  if (jSurfaces.empty() && jVolumes.empty()) {
    auto volume = Experimental::DetectorVolumeFactory::construct(
        portalGenerator, gctx, name, transform, std::move(bounds),
        Experimental::tryAllPortals());
    volume->assignGeometryId(geoId);
    return volume;
  }
  // Convert the surfaces
  std::vector<std::shared_ptr<Surface>> surfaces;
  for (const auto& js : jSurfaces) {
    surfaces.push_back(SurfaceJsonConverter::fromJson(js));
  }
  // Convert the volumes
  std::vector<std::shared_ptr<Experimental::DetectorVolume>> volumes;
  for (const auto& jv : jVolumes) {
    volumes.push_back(DetectorVolumeJsonConverter::fromJson(gctx, jv));
  }

  auto jSurfaceNavigation = jVolume["surface_navigation"];

  auto volume = Experimental::DetectorVolumeFactory::construct(
      portalGenerator, gctx, name, transform, std::move(bounds), surfaces,
      volumes, Experimental::tryRootVolumes(),
      IndexedSurfacesJsonConverter::fromJson(jSurfaceNavigation));
  volume->assignGeometryId(geoId);
  return volume;
}
