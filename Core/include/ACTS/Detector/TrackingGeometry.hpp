// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingGeometry.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_TRACKINGGEOMETRY_H
#define ACTS_DETECTOR_TRACKINGGEOMETRY_H

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometrySignature.hpp"
// STD
#include <map>
#include <memory>
#include <vector>

namespace Acts {

class TrackingVolume;
class DetachedTrackingVolume;
class PerigeeSurface;
class Layer;

typedef std::shared_ptr<const TrackingVolume>         TrackingVolumePtr;
typedef std::shared_ptr<const DetachedTrackingVolume> DetachedTrackingVolumePtr;
typedef std::vector<DetachedTrackingVolumePtr>        DetachedVolumeVector;

///  @class TrackingGeometry
///
///  The TrackingGeometry class is the owner of the constructed TrackingVolumes.
///
///  It enables both, a global search for an asociatedVolume
///  (respectively, if existing, a global search of an associated Layer or the next
///  associated Layer), such as a continous navigation by BoundarySurfaces between
///  the confined TrackingVolumes.

class TrackingGeometry
{
  /// Give the GeometryBuilder friend rights 
  friend class TrackingGeometryBuilder;

public:
  /// Constructor
  TrackingGeometry(TrackingVolumePtr highestVolume);

  /// Destructor 
  ~TrackingGeometry();

  /// return the world 
  const TrackingVolume*
  highestTrackingVolume() const;

  /// return the lowest tracking Volume
  /// @param gpos is the global position fo the call
  const TrackingVolume*
  lowestTrackingVolume(const Vector3D& gpos) const;

  /// return the vector of lowest detached tracking Volume(->overlaps)
  /// @param gpos is the global position fo the call
  const DetachedVolumeVector*
  lowestDetachedTrackingVolumes(const Vector3D& gpos) const;

  /// return the lowest static volume
  /// @param gpos is the global position fo the call
  const TrackingVolume*
  lowestStaticTrackingVolume(const Vector3D& gpos) const;

  /// return the tracking Volume by name, 0 if it doesn't exist 
  const TrackingVolume*
  trackingVolume(const std::string& name) const;

  /// Forward the associated Layer information 
  /// @param gpos is the global position fo the call
  const Layer*
  associatedLayer(const Vector3D& gpos) const;

  /// check position at volume boundary 
  /// @param gpos is the global position fo the call
  bool
  atVolumeBoundary(const Vector3D&       gpos,
                   const TrackingVolume* vol,
                   double                tol) const;

  /// check position at volume boundary + navigation link 
  bool
  atVolumeBoundary(const Vector3D&        gpos,
                   const Vector3D&        mom,
                   const TrackingVolume*  vol,
                   const TrackingVolume*& nextVol,
                   Acts::PropDirection    dir,
                   double                 tol) const;
 
  /// register the beam tube 
  void
  registerBeamTube(std::unique_ptr<const PerigeeSurface> beam) const;

private:
  /// Geometry Builder busineess: the geometry builder has to sign
  void
  sign(GeometrySignature geosit, GeometryType geotype = Static) const;

  /// private method to register recursively the tracking volume
  /// & set the mother volume 
  void
  registerTrackingVolumes(const TrackingVolume& tvol,
                          const TrackingVolume* mvol = nullptr,
                          int                   lvl  = 0);

  /// The known world - and the beamline 
  TrackingVolumePtr                             m_world;
  mutable std::unique_ptr<const PerigeeSurface> m_beam;

  /// The Volumes in a map for later finding 
  std::map<const std::string, const TrackingVolume*> m_trackingVolumes;
};

}  // end of namespace

#endif
