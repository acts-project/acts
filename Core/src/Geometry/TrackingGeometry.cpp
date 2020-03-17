// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingGeometry.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include <functional>

#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

Acts::TrackingGeometry::TrackingGeometry(
    const MutableTrackingVolumePtr& highestVolume,
    const IMaterialDecorator* materialDecorator)
    : m_world(highestVolume),
      m_beam(Surface::makeShared<PerigeeSurface>(s_origin)) {
  // Close the geometry: assign geometryID and successively the material
  size_t volumeID = 0;
  highestVolume->closeGeometry(materialDecorator, m_trackingVolumes, volumeID);
}

Acts::TrackingGeometry::~TrackingGeometry() = default;

const Acts::TrackingVolume* Acts::TrackingGeometry::lowestTrackingVolume(
    const GeometryContext& gctx, const Acts::Vector3D& gp) const {
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
  return (m_world.get());
}

const Acts::TrackingVolume* Acts::TrackingGeometry::trackingVolume(
    const std::string& name) const {
  auto sVol = m_trackingVolumes.begin();
  sVol = m_trackingVolumes.find(name);
  if (sVol != m_trackingVolumes.end()) {
    return (sVol->second);
  }
  return nullptr;
}

const Acts::Layer* Acts::TrackingGeometry::associatedLayer(
    const GeometryContext& gctx, const Acts::Vector3D& gp) const {
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

void Acts::TrackingGeometry::visitSurfaces(
    const std::function<void(const Acts::Surface*)>& visitor) const {
  highestTrackingVolume()->visitSurfaces(visitor);
}
