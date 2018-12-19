// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingGeometry.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometrySignature.hpp"
// STD
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace Acts {

class TrackingVolume;
class DetachedTrackingVolume;
class Surface;
class PerigeeSurface;
class Layer;

using TrackingVolumePtr         = std::shared_ptr<const TrackingVolume>;
using MutableTrackingVolumePtr  = std::shared_ptr<TrackingVolume>;
using DetachedTrackingVolumePtr = std::shared_ptr<const DetachedTrackingVolume>;
using DetachedVolumeVector      = std::vector<DetachedTrackingVolumePtr>;

///  @class TrackingGeometry
///
///  The TrackingGeometry class is the owner of the constructed TrackingVolumes.
///
///  It enables both, a global search for an asociatedVolume
///  (respectively, if existing, a global search of an associated Layer or the
///  next associated Layer), such as a continous navigation by BoundarySurfaces
///  between the confined TrackingVolumes.

class TrackingGeometry
{
  /// Give the GeometryBuilder friend rights
  friend class TrackingGeometryBuilder;

public:
  /// Constructor
  ///
  /// @param highestVolume is the world volume
  TrackingGeometry(const MutableTrackingVolumePtr& highestVolume);

  /// Destructor
  ~TrackingGeometry();

  /// Access to the world volume
  /// @return plain pointer to the world volume
  const TrackingVolume*
  highestTrackingVolume() const;

  /// return the lowest tracking Volume
  ///
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to the lowest TrackingVolume
  const TrackingVolume*
  lowestTrackingVolume(const Vector3D& gp) const;

  /// return the vector of lowest detached tracking Volume(->overlaps)
  ///
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to the the lowest DetachedTrackingVolume
  const DetachedVolumeVector*
  lowestDetachedTrackingVolumes(const Vector3D& gp) const;

  /// return the lowest static volume
  ///
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to the the lowest static tracking volume
  const TrackingVolume*
  lowestStaticTrackingVolume(const Vector3D& gp) const;

  /// return the lowest tracking Volume
  ///
  /// @param name is the name for the volume search
  ///
  /// @return plain pointer to the lowest TrackingVolume
  const TrackingVolume*
  trackingVolume(const std::string& name) const;

  /// Forward the associated Layer information
  ///
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to assocaiated layer
  const Layer*
  associatedLayer(const Vector3D& gp) const;

  /// Register the beam tube
  ///
  /// @param beam is the beam line surface
  void
  registerBeamTube(std::shared_ptr<const PerigeeSurface> beam);

  /// @brief surface representing the beam pipe
  ///
  /// @note The ownership is not passed, e.g. do not delete the pointer
  ///
  /// @return raw pointer to surface representing the beam pipe
  ///         (could be a null pointer)
  const Surface*
  getBeamline() const;

  /// @brief Visit all sensitive surfaces
  ///
  /// @param visitor The callable. Will be called for each sensitive surface
  /// that is found
  void
  visitSurfaces(const std::function<void(const Acts::Surface*)>& visitor) const;

private:
  /// Geometry Builder busineess: the geometry builder has to sign
  ///
  /// @param geosit is the volume signature
  /// @param geotype is the volume navigation type
  void
  sign(GeometrySignature geosit, GeometryType geotype = Static);

  /// The known world - and the beamline
  TrackingVolumePtr                     m_world;
  std::shared_ptr<const PerigeeSurface> m_beam;

  /// The Volumes in a map for string based search
  std::map<std::string, const TrackingVolume*> m_trackingVolumes;
};

}  // namespace
