// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingGeometry.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/DetachedTrackingVolume.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/PerigeeSurface.hpp"
#include "ACTS/Surfaces/Surface.hpp"

Acts::TrackingGeometry::TrackingGeometry(TrackingVolumePtr highestVolume)
  : m_world(highestVolume)
  , m_beam(std::make_unique<const PerigeeSurface>(s_origin))
{
  // create the GeometryID for this
  GeometryID geoID(0);
  // close up the geometry
  if (m_world) m_world->closeGeometry(geoID, m_trackingVolumes);
}

Acts::TrackingGeometry::~TrackingGeometry()
{
}

const Acts::TrackingVolume*
Acts::TrackingGeometry::lowestTrackingVolume(const Acts::Vector3D& gp) const
{
  const TrackingVolume* searchVolume  = m_world.get();
  const TrackingVolume* currentVolume = nullptr;
  while (currentVolume != searchVolume && searchVolume) {
    currentVolume = searchVolume;
    searchVolume  = searchVolume->trackingVolume(gp);
  }
  return currentVolume;
}

const Acts::DetachedVolumeVector*
Acts::TrackingGeometry::lowestDetachedTrackingVolumes(const Vector3D& gp) const
{
  double                      tol           = 0.001;
  const TrackingVolume*       currentVolume = lowestStaticTrackingVolume(gp);
  if (currentVolume) return currentVolume->detachedTrackingVolumes(gp, tol);
  return nullptr;
}

const Acts::TrackingVolume*
Acts::TrackingGeometry::lowestStaticTrackingVolume(const Vector3D& gp) const
{
  const TrackingVolume* searchVolume  = m_world.get();
  const TrackingVolume* currentVolume = nullptr;
  while (currentVolume != searchVolume && searchVolume) {
    currentVolume = searchVolume;
    if ((searchVolume->confinedDetachedVolumes()).empty())
      searchVolume = searchVolume->trackingVolume(gp);
  }
  return currentVolume;
}

//@todo change to BoundaryCheck
bool
Acts::TrackingGeometry::atVolumeBoundary(const Acts::Vector3D&       gp,
                                         const Acts::TrackingVolume* vol,
                                         double) const
{
  bool isAtBoundary = false;
  if (!vol) return isAtBoundary;
  for (auto& bSurface : vol->boundarySurfaces()) {
    const Surface& surf = bSurface->surfaceRepresentation();
    if (surf.isOnSurface(gp, true)) isAtBoundary = true;
  }
  return isAtBoundary;
}

//@todo change to BoundaryCheck
bool
Acts::TrackingGeometry::atVolumeBoundary(const Vector3D&        gp,
                                         const Vector3D&        mom,
                                         const TrackingVolume*  vol,
                                         const TrackingVolume*& nextVol,
                                         PropDirection          dir,
                                         double) const
{
  bool isAtBoundary = false;
  nextVol           = 0;
  if (!vol) return isAtBoundary;
  for (auto& bSurface : vol->boundarySurfaces()) {
    const Surface& surf = bSurface->surfaceRepresentation();
    if (surf.isOnSurface(gp, true)) {
      isAtBoundary = true;
      const TrackingVolume* attachedVol
          = bSurface->attachedVolume(gp, mom, dir);
      if (!nextVol && attachedVol) nextVol = attachedVol;
    }
  }
  return isAtBoundary;
}

const Acts::TrackingVolume*
Acts::TrackingGeometry::highestTrackingVolume() const
{
  return (m_world.get());
}

void
Acts::TrackingGeometry::sign(GeometrySignature geosit,
                             GeometryType      geotype)
{
  auto mutableWorld = std::const_pointer_cast<TrackingVolume>( m_world );
  mutableWorld->sign(geosit, geotype);
}

const Acts::TrackingVolume*
Acts::TrackingGeometry::trackingVolume(const std::string& name) const
{
  auto sVol = m_trackingVolumes.begin();
  sVol      = m_trackingVolumes.find(name);
  if (sVol != m_trackingVolumes.end()) return (sVol->second);
  return nullptr;
}

const Acts::Layer*
Acts::TrackingGeometry::associatedLayer(const Acts::Vector3D& gp) const
{
  const TrackingVolume* lowestVol = (lowestTrackingVolume(gp));
  return lowestVol->associatedLayer(gp);
}

void
Acts::TrackingGeometry::registerBeamTube(
    std::unique_ptr<const PerigeeSurface> beam)
{
  m_beam = std::move(beam);
}

const Acts::Surface*
Acts::TrackingGeometry::getBeamline() const
{
  return m_beam.get();
}
