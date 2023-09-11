// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/SequentialGeometryIdGenerator.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"

Acts::Experimental::IGeometryIdGenerator::GeoIdCache
Acts::Experimental::SequentialGeometryIdGenerator::generateCache() const {
  return Cache{};
}

void Acts::Experimental::SequentialGeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, DetectorVolume& dVolume) const {
  auto ccache = std::any_cast<Cache&>(cache);

  // Set to the volume itself
  if (dVolume.geometryId().volume() == 0 or m_cfg.overrideExistingIds) {
    GeometryIdentifier geoID(0u);
    geoID.setVolume(++ccache.volumeCount);
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

void Acts::Experimental::SequentialGeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, Portal& portal) const {
  auto ccache = std::any_cast<Cache&>(cache);

  auto& pSurface = portal.surface();
  if (pSurface.geometryId().boundary() == 0 or m_cfg.overrideExistingIds) {
    GeometryIdentifier geoID(0u);
    geoID.setVolume(ccache.volumeCount);
    geoID.setBoundary(++ccache.portalCount);
    pSurface.assignGeometryId(geoID);
  }
}

void Acts::Experimental::SequentialGeometryIdGenerator::assignGeometryId(
    IGeometryIdGenerator::GeoIdCache& cache, Surface& surface) const {
  auto ccache = std::any_cast<Cache&>(cache);

  auto rGeoID = surface.geometryId();
  if ((rGeoID.sensitive() == 0 and rGeoID.passive() == 0) or
      m_cfg.overrideExistingIds) {
    GeometryIdentifier geoID(0u);
    geoID.setVolume(ccache.volumeCount);
    if (surface.associatedDetectorElement() != nullptr) {
      geoID.setSensitive(++ccache.sensitiveCount);
    } else {
      geoID.setPassive(++ccache.passiveCount);
    }
    surface.assignGeometryId(geoID);
  }
}
