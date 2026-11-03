// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/detail/MaterialSurfaceRegistry.hpp"

#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <functional>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace Acts::detail {

const std::optional<std::string>& materialKey(
    const ISurfaceMaterial* material) {
  if (auto* proto = dynamic_cast<const ProtoSurfaceMaterial*>(material)) {
    return proto->materialKey();
  }
  if (auto* proto = dynamic_cast<const ProtoGridSurfaceMaterial*>(material)) {
    return proto->materialKey();
  }
  static const std::optional<std::string> empty;
  return empty;
}

MaterialSurfaceRegistry::MaterialSurfaceRegistry(
    std::span<const Surface* const> surfaces) {
  std::map<std::string, const Surface*, std::less<>> keys;
  std::unordered_set<const Surface*> visited;
  for (const Surface* surface : surfaces) {
    if (surface == nullptr) {
      throw std::invalid_argument("Null material surface");
    }
    if (!visited.insert(surface).second) {
      continue;
    }
    if (auto* marker = dynamic_cast<const MergedMaterialMarker*>(
            surface->surfaceMaterial())) {
      std::ostringstream message;
      message << "Cannot map or load material on surface "
              << surface->geometryId() << ": " << *marker;
      throw std::invalid_argument(message.str());
    }
    const auto id = surface->geometryId();
    if (!m_surfaces.emplace(id, surface).second) {
      std::ostringstream message;
      message << "Duplicate material surface geometry ID " << id;
      throw std::invalid_argument(message.str());
    }
    m_keys.try_emplace(id, detail::materialKey(surface->surfaceMaterial()));
    if (detail::materialKey(surface->surfaceMaterial())) {
      const auto& key = *detail::materialKey(surface->surfaceMaterial());
      if (id == GeometryIdentifier{}) {
        throw std::invalid_argument("Material key '" + key +
                                    "' has no geometry ID");
      }
      if (!keys.try_emplace(key, surface).second) {
        throw std::invalid_argument("Duplicate material key '" + key + "'");
      }
    }
  }
}

TrackingGeometryMaterial MaterialSurfaceRegistry::materialMaps(
    SurfaceMaterialMaps materials) const {
  TrackingGeometryMaterial result;
  for (const auto& [id, key] : m_keys) {
    if (key && !materials.contains(id)) {
      throw std::invalid_argument("Missing mapping output for material key '" +
                                  *key + "'");
    }
  }
  for (auto& [id, material] : materials) {
    const auto found = m_surfaces.find(id);
    if (found == m_surfaces.end()) {
      throw std::invalid_argument(
          "Material output has an unregistered geometry ID");
    }
    const Surface& surface = *found->second;
    if (surface.geometryId() != id ||
        detail::materialKey(surface.surfaceMaterial()) != m_keys.at(id)) {
      throw std::invalid_argument(
          "Material surface identity changed during mapping");
    }
    if (!material) {
      throw std::invalid_argument("Null material in mapping output");
    }
    if (detail::materialKey(surface.surfaceMaterial())) {
      if (!result.keyedSurfaces
               .try_emplace(*detail::materialKey(surface.surfaceMaterial()),
                            KeyedSurfaceMaterial{id, std::move(material)})
               .second) {
        throw std::invalid_argument(
            "Duplicate material key during finalization");
      }
    } else {
      result.surfaceMaterials.try_emplace(id, std::move(material));
    }
  }
  return result;
}

}  // namespace Acts::detail
