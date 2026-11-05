// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"

#include <format>
#include <stdexcept>

namespace ActsPlugins {

DetraySurfaceIdMap::DetraySurfaceIdMap(std::span<const Entry> entries) {
  m_fromDetray.reserve(entries.size());
  m_toDetray.reserve(entries.size());
  for (const auto& [detrayId, surface] : entries) {
    if (surface == nullptr) {
      throw std::invalid_argument(
          std::format("DetraySurfaceIdMap: null surface for detray id {}",
                      detrayId.value()));
    }
    const Acts::GeometryIdentifier id = surface->geometryId();
    if (m_fromDetray.contains(detrayId) || m_toDetray.contains(id)) {
      throw std::invalid_argument(
          std::format("DetraySurfaceIdMap: duplicate surface (detray id {}, "
                      "ACTS id {})",
                      detrayId.value(), id));
    }
    m_fromDetray.emplace(detrayId, surface);
    m_toDetray.emplace(id, detrayId);
  }
}

std::optional<detray::geometry::identifier> DetraySurfaceIdMap::toDetray(
    Acts::GeometryIdentifier id) const {
  auto it = m_toDetray.find(id);
  if (it == m_toDetray.end()) {
    return std::nullopt;
  }
  return it->second;
}

std::optional<Acts::GeometryIdentifier> DetraySurfaceIdMap::fromDetray(
    detray::geometry::identifier detrayId) const {
  const Acts::Surface* actsSurface = surface(detrayId);
  if (actsSurface == nullptr) {
    return std::nullopt;
  }
  return actsSurface->geometryId();
}

const Acts::Surface* DetraySurfaceIdMap::surface(
    detray::geometry::identifier detrayId) const {
  auto it = m_fromDetray.find(detrayId);
  return it == m_fromDetray.end() ? nullptr : it->second;
}

}  // namespace ActsPlugins
