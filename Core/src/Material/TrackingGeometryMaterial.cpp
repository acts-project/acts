// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/TrackingGeometryMaterial.hpp"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/detail/MaterialSurfaceRegistry.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Acts {

namespace {
std::shared_ptr<const ISurfaceMaterial> resolve(
    const TrackingGeometryMaterial& maps, const Surface& surface) {
  if (const auto* marker = dynamic_cast<const MergedMaterialMarker*>(
          surface.surfaceMaterial())) {
    std::ostringstream message;
    message << "Cannot decorate surface " << surface.geometryId() << ": "
            << *marker;
    throw std::invalid_argument(message.str());
  }
  if (detail::materialKey(surface.surfaceMaterial())) {
    const auto& key = *detail::materialKey(surface.surfaceMaterial());
    const auto found = maps.keyedSurfaces.find(key);
    if (found == maps.keyedSurfaces.end() || !found->second.material) {
      throw std::invalid_argument("Missing material for key '" + key + "'");
    }
    auto checkProtoKey = [&](const auto* proto) {
      if (proto != nullptr && proto->materialKey() &&
          *proto->materialKey() != key) {
        throw std::invalid_argument(
            "Proto material disagrees with assignment key '" + key + "'");
      }
    };
    checkProtoKey(dynamic_cast<const ProtoSurfaceMaterial*>(
        found->second.material.get()));
    checkProtoKey(dynamic_cast<const ProtoGridSurfaceMaterial*>(
        found->second.material.get()));
    if (dynamic_cast<const MergedMaterialMarker*>(
            found->second.material.get()) != nullptr) {
      throw std::invalid_argument("Merged material marker cannot supply key '" +
                                  key + "'");
    }
    return found->second.material;
  }
  const auto found = maps.surfaceMaterials.find(surface.geometryId());
  return found == maps.surfaceMaterials.end() ? nullptr : found->second;
}

}  // namespace

void TrackingGeometryMaterial::apply(Surface& surface) const {
  if (auto material = resolve(*this, surface)) {
    surface.assignSurfaceMaterial(std::move(material));
  }
}

void TrackingGeometryMaterial::apply(std::span<Surface* const> surfaces) const {
  std::vector<const Surface*> participating;
  for (const Surface* surface : surfaces) {
    if (surface == nullptr) {
      throw std::invalid_argument("Null surface during material decoration");
    }
    if (surface->hasMaterial() ||
        detail::materialKey(surface->surfaceMaterial()) ||
        surfaceMaterials.contains(surface->geometryId())) {
      participating.push_back(surface);
    }
  }
  detail::MaterialSurfaceRegistry registry(participating);
  std::vector<std::pair<Surface*, std::shared_ptr<const ISurfaceMaterial>>>
      assignments;
  std::set<Surface*> visited;
  for (Surface* surface : surfaces) {
    if (visited.insert(surface).second) {
      if (auto material = resolve(*this, *surface)) {
        assignments.emplace_back(surface, std::move(material));
      }
    }
  }
  for (auto& [surface, material] : assignments) {
    surface->assignSurfaceMaterial(std::move(material));
  }
}

void TrackingGeometryMaterial::apply(TrackingGeometry& geometry) const {
  if (!keyedSurfaces.empty() &&
      geometry.geometryVersion() == TrackingGeometry::GeometryVersion::Gen1) {
    throw std::invalid_argument(
        "Cannot apply a keyed material map to Gen1 geometry: stable material "
        "keys require Gen3 material designators. Use a Gen1 map indexed by "
        "geometry ID or apply this map to the matching Gen3 geometry.");
  }
  std::vector<Surface*> surfaces;
  std::vector<TrackingVolume*> volumes;
  geometry.apply(
      overloaded{[&](Surface& surface) { surfaces.push_back(&surface); },
                 [&](Portal& portal) { surfaces.push_back(&portal.surface()); },
                 [&](BoundarySurfaceT<TrackingVolume>& boundary) {
                   surfaces.push_back(&boundary.surfaceRepresentation());
                 },
                 [&](TrackingVolume& volume) { volumes.push_back(&volume); }});
  apply(surfaces);
  for (auto* volume : volumes) {
    apply(*volume);
  }
}

void TrackingGeometryMaterial::apply(TrackingVolume& volume) const {
  if (const auto found = volumeMaterials.find(volume.geometryId());
      found != volumeMaterials.end()) {
    volume.assignVolumeMaterial(found->second);
  }
}

}  // namespace Acts
