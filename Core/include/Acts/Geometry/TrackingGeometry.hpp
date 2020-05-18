// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace Acts {

class TrackingVolume;
class Layer;
class Surface;
class PerigeeSurface;
class IMaterialDecorator;

using TrackingVolumePtr = std::shared_ptr<const TrackingVolume>;
using MutableTrackingVolumePtr = std::shared_ptr<TrackingVolume>;

///  @class TrackingGeometry
///
///  The TrackingGeometry class is the owner of the constructed TrackingVolumes.
///
///  It enables both, a global search for an asociatedVolume
///  (respectively, if existing, a global search of an associated Layer or the
///  next associated Layer), such as a continous navigation by BoundarySurfaces
///  between the confined TrackingVolumes.
class TrackingGeometry {
  /// Give the GeometryBuilder friend rights
  friend class TrackingGeometryBuilder;

 public:
  /// Constructor
  ///
  /// @param highestVolume is the world volume
  /// @param materialDecorator is a dediated decorator that can assign
  ///        surface or volume based material to the TrackingVolume
  TrackingGeometry(const MutableTrackingVolumePtr& highestVolume,
                   const IMaterialDecorator* materialDecorator = nullptr);

  /// Destructor
  ~TrackingGeometry();

  /// Access to the world volume
  /// @return plain pointer to the world volume
  const TrackingVolume* highestTrackingVolume() const;

  /// return the lowest tracking Volume
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to the lowest TrackingVolume
  const TrackingVolume* lowestTrackingVolume(const GeometryContext& gctx,
                                             const Vector3D& gp) const;

  /// return the lowest tracking Volume
  ///
  /// @param name is the name for the volume search
  ///
  /// @return plain pointer to the lowest TrackingVolume
  const TrackingVolume* trackingVolume(const std::string& name) const;

  /// Forward the associated Layer information
  ///
  /// @paramn gctx is the context for this request (e.g. alignment)
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to assocaiated layer
  const Layer* associatedLayer(const GeometryContext& gctx,
                               const Vector3D& gp) const;

  /// Register the beam tube
  ///
  /// @param beam is the beam line surface
  void registerBeamTube(std::shared_ptr<const PerigeeSurface> beam);

  /// @brief surface representing the beam pipe
  ///
  /// @note The ownership is not passed, e.g. do not delete the pointer
  ///
  /// @return raw pointer to surface representing the beam pipe
  ///         (could be a null pointer)
  const Surface* getBeamline() const;

  /// @brief Visit all sensitive surfaces
  ///
  /// @param visitor The callable. Will be called for each sensitive surface
  /// that is found
  void visitSurfaces(
      const std::function<void(const Acts::Surface*)>& visitor) const;

 private:
  /// The known world - and the beamline
  TrackingVolumePtr m_world;
  std::shared_ptr<const PerigeeSurface> m_beam;

  /// The Volumes in a map for string based search
  std::map<std::string, const TrackingVolume*> m_trackingVolumes;
};

}  // namespace Acts
