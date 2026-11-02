// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"

#include <sstream>
#include <stdexcept>

namespace ActsPlugins {

void DetraySurfaceIdMap::insert(DetrayId detrayId,
                                Acts::GeometryIdentifier id) {
  if (m_fromDetray.contains(detrayId) || m_toDetray.contains(id)) {
    std::stringstream msg;
    msg << "DetraySurfaceIdMap: duplicate surface (detray id " << detrayId
        << ", ACTS id " << id << ")";
    throw std::invalid_argument(msg.str());
  }
  m_fromDetray.emplace(detrayId, id);
  m_toDetray.emplace(id, detrayId);
}

std::optional<DetraySurfaceIdMap::DetrayId> DetraySurfaceIdMap::toDetray(
    Acts::GeometryIdentifier id) const {
  auto it = m_toDetray.find(id);
  if (it == m_toDetray.end()) {
    return std::nullopt;
  }
  return it->second;
}

std::optional<Acts::GeometryIdentifier> DetraySurfaceIdMap::fromDetray(
    DetrayId detrayId) const {
  auto it = m_fromDetray.find(detrayId);
  if (it == m_fromDetray.end()) {
    return std::nullopt;
  }
  return it->second;
}

}  // namespace ActsPlugins
