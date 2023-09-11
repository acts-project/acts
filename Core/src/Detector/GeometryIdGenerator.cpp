// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/GeometryIdGenerator.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Surfaces/Surface.hpp"

Acts::Experimental::IGeometryIdGenerator::GeoIdCache
Acts::Experimental::GeometryIdGenerator::generateCache() const {
  return Cache{};
}

void Acts::Experimental::GeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, DetectorVolume& dVolume) const {
  auto ccache = std::any_cast<Cache&>(cache);

  // Set to the volume itself
  if (dVolume.geometryId().volume() == 0 or m_cfg.overrideExistingIds) {
    ++ccache.volumeCount;
    GeometryIdentifier geoID = volumeId(ccache);
    dVolume.assignGeometryId(geoID);
  }

  // Portals
  std::for_each(dVolume.portalPtrs().begin(), dVolume.portalPtrs().end(),
                [&](auto& portal) { assignGeometryId(cache, *portal); });

  // Surfaces
  std::for_each(dVolume.surfacePtrs().begin(), dVolume.surfacePtrs().end(),
                [&](auto& surface) { assignGeometryId(cache, *surface); });

  // Sub volumes
  std::for_each(dVolume.volumePtrs().begin(), dVolume.volumePtrs().end(),
                [&](auto& volume) { assignGeometryId(cache, *volume); });
}

void Acts::Experimental::GeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, Portal& portal) const {
  auto ccache = std::any_cast<Cache&>(cache);

  auto& pSurface = portal.surface();
  if (pSurface.geometryId().boundary() == 0 or m_cfg.overrideExistingIds) {
    GeometryIdentifier geoID = volumeId(ccache);
    geoID.setBoundary(++ccache.portalCount);
    pSurface.assignGeometryId(geoID);
  }
}

void Acts::Experimental::GeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, Surface& surface) const {
  auto ccache = std::any_cast<Cache&>(cache);

  auto rGeoID = surface.geometryId();
  if ((rGeoID.sensitive() == 0 and rGeoID.passive() == 0) or
      m_cfg.overrideExistingIds) {
    GeometryIdentifier geoID = volumeId(ccache);
    if (surface.associatedDetectorElement() != nullptr) {
      geoID.setSensitive(++ccache.sensitiveCount);
    } else {
      geoID.setPassive(++ccache.passiveCount);
    }
    surface.assignGeometryId(geoID);
  }
}

Acts::GeometryIdentifier Acts::Experimental::GeometryIdGenerator::volumeId(
    Cache& cache) const {
  GeometryIdentifier geoID(0u);
  if (not m_cfg.containerMode) {
    geoID.setVolume(cache.volumeCount);
  } else {
    geoID.setVolume(m_cfg.containerId);
    geoID.setLayer(++cache.layerCount);
  }
  return geoID;
}
