// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

Acts::TrackingGeometry::TrackingGeometry(
    const MutableTrackingVolumePtr& highestVolume,
    const IMaterialDecorator* materialDecorator,
    const GeometryIdentifierHook& hook, const Logger& logger)
    : m_world(highestVolume),
      m_beam(Surface::makeShared<PerigeeSurface>(Vector3::Zero())) {
  // Close the geometry: assign geometryID and successively the material
  std::size_t volumeID = 0;
  highestVolume->closeGeometry(materialDecorator, m_volumesById, volumeID, hook,
                               logger);
  m_volumesById.rehash(0);
  // fill surface lookup container
  m_world->visitSurfaces([this](const Acts::Surface* srf) {
    if (srf != nullptr) {
      m_surfacesById[srf->geometryId()] = srf;
    }
  });
  m_surfacesById.rehash(0);
}

Acts::TrackingGeometry::~TrackingGeometry() = default;

const Acts::TrackingVolume* Acts::TrackingGeometry::lowestTrackingVolume(
    const GeometryContext& gctx, const Acts::Vector3& gp) const {
  const TrackingVolume* searchVolume = m_world.get();
  const TrackingVolume* currentVolume = nullptr;
  while (currentVolume != searchVolume && (searchVolume != nullptr)) {
    currentVolume = searchVolume;
    searchVolume = searchVolume->lowestTrackingVolume(gctx, gp);
  }
  return currentVolume;
}

const Acts::TrackingVolume* Acts::TrackingGeometry::highestTrackingVolume()
    const {
  return m_world.get();
}

const std::shared_ptr<const Acts::TrackingVolume>&
Acts::TrackingGeometry::highestTrackingVolumeShared() const {
  return m_world;
}

const Acts::Layer* Acts::TrackingGeometry::associatedLayer(
    const GeometryContext& gctx, const Acts::Vector3& gp) const {
  const TrackingVolume* lowestVol = (lowestTrackingVolume(gctx, gp));
  return lowestVol->associatedLayer(gctx, gp);
}

void Acts::TrackingGeometry::registerBeamTube(
    std::shared_ptr<const PerigeeSurface> beam) {
  m_beam = std::move(beam);
}

const Acts::Surface* Acts::TrackingGeometry::getBeamline() const {
  return m_beam.get();
}

const Acts::TrackingVolume* Acts::TrackingGeometry::findVolume(
    GeometryIdentifier id) const {
  auto vol = m_volumesById.find(id);
  if (vol == m_volumesById.end()) {
    return nullptr;
  }
  return vol->second;
}

const Acts::Surface* Acts::TrackingGeometry::findSurface(
    GeometryIdentifier id) const {
  auto srf = m_surfacesById.find(id);
  if (srf == m_surfacesById.end()) {
    return nullptr;
  }
  return srf->second;
}

const std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>&
Acts::TrackingGeometry::geoIdSurfaceMap() const {
  return m_surfacesById;
}
