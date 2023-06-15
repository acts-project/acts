// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/SurfaceVisitorConcept.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace Acts {

class Layer;
class Surface;
class PerigeeSurface;
class IMaterialDecorator;
class TrackingVolume;

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
  /// @param hook Identifier hook to be applied to surfaces
  /// @param logger instance of a logger (defaulting to the "silent" one)
  TrackingGeometry(const MutableTrackingVolumePtr& highestVolume,
                   const IMaterialDecorator* materialDecorator = nullptr,
                   const GeometryIdentifierHook& hook = {},
                   const Logger& logger = getDummyLogger());

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
                                             const Vector3& gp) const;

  /// Forward the associated Layer information
  ///
  /// @param gctx is the context for this request (e.g. alignment)
  /// @param gp is the global position of the call
  ///
  /// @return plain pointer to assocaiated layer
  const Layer* associatedLayer(const GeometryContext& gctx,
                               const Vector3& gp) const;

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
  /// @tparam visitor_t Type of the callable visitor
  ///
  /// @param visitor The callable. Will be called for each sensitive surface
  /// that is found
  template <ACTS_CONCEPT(SurfaceVisitor) visitor_t>
  void visitSurfaces(visitor_t&& visitor) const {
    highestTrackingVolume()->template visitSurfaces<visitor_t>(
        std::forward<visitor_t>(visitor));
  }

  /// Search for a volume with the given identifier.
  ///
  /// @param id is the geometry identifier of the volume
  /// @retval nullptr if no such volume exists
  /// @retval pointer to the found volume otherwise.
  const TrackingVolume* findVolume(GeometryIdentifier id) const;

  /// Search for a surface with the given identifier.
  ///
  /// @param id is the geometry identifier of the surface
  /// @retval nullptr if no such surface exists
  /// @retval pointer to the found surface otherwise.
  const Surface* findSurface(GeometryIdentifier id) const;

 private:
  // the known world
  TrackingVolumePtr m_world;
  // beam line
  std::shared_ptr<const PerigeeSurface> m_beam;
  // lookup containers
  std::unordered_map<GeometryIdentifier, const TrackingVolume*> m_volumesById;
  std::unordered_map<GeometryIdentifier, const Surface*> m_surfacesById;
};

}  // namespace Acts
