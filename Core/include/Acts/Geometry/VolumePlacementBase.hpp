// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/detail/PortalPlacement.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"

#include <memory>
#include <vector>

namespace Acts {
class GeometryContext;

/// @brief Interface class to define the backend for alignable volumes
///        that move coherently with the snesitive surfaces inside
///        The interface provides the transform from local -> global
///        coordinates as well as its way back.
///        To move the oriented surfaces along with the volume itself
///        the interface needs to provide a factory that creates the
///        detector elements that are then passed to the oriented surfaces
class VolumePlacementBase {
 public:
  /// @brief Default constructor
  VolumePlacementBase();
  /// @brief Default destructor
  virtual ~VolumePlacementBase() = default;
  /// @brief Attaches a SurfacePlacementBase to the volume's portal surface
  ///        to synchronize the portal alignment with the volume alignment
  ///        Internally, a `PortalPlacement` object is instantiated and
  ///        connected with the surface. As the number of calls to create
  ///        the portal surfaces is a priori unrestricted, the ownership
  ///        of the portal surface is transferred
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  /// @param portalToVolTrf: Transform to switch from the portal's frame into the volume's frame
  /// @param surface: Mutable reference to the portal surface that is to be aligned
  std::shared_ptr<RegularSurface> makePortalAlignable(
      const std::size_t portalIdx, const Transform3& portalToVolTrf,
      std::shared_ptr<RegularSurface>&& surface);
  /// @brief Returns the const pointer to the `SurfacePlacementBase` object
  ///        aligning the i-th portal (May be nullptr if index exceeds the
  ///        number of portals)
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  const detail::PortalPlacement* portalPlacement(
      const std::size_t portalIdx) const;
  /// @brief Returns the number of portal placement objects
  std::size_t nPortalPlacement() const;
  /// @brief Returns the transformation from the local volume coordinates to
  ///        the experiment's global coordinate system
  /// @param gctx The current geometry context object, e.g. alignment
  virtual const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const = 0;
  /// @brief Returns the transformation from the experiment's global frame to the
  ///        local volume coordinate system
  /// @param gctx The current geometry context object, e.g. alignment
  virtual const Transform3& globalToLocalTransform(
      const GeometryContext& gctx) const = 0;
  /// @brief Returns the transform from the portal's frame to the experiment's
  ///        global frame for the portal surface associated with the volume
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  virtual const Transform3& portalLocalToGlobal(
      const GeometryContext& gctx, const std::size_t portalIdx) const = 0;
  /// @brief Helper method to invoke the population of the experiment specific
  ///        GeometryContext with the transforms of the associated portals
  ///        If the implementation of the GeometryContext does not forsee
  ///        a const-mutable lazy population, this method has to be called
  ///        at every change of the alignment constants
  /// @param gctx: Mutable reference of the GeometryContext to populate
  void populateContextWithPortals(GeometryContext& gctx) const;

 protected:
  /// @brief Invoke the backend cache to store the localToGlobalTransform
  ///        for the given portal in the geometry Context
  /// @param gctx: Reference to the mutable GeometryContext where the
  ///              tarnsform for the portal shall be stored
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  /// @param portalLocToGlob: Transform to switch from the portal's frame to the
  ///                         experiment's global coordinate frame
  virtual void cachePortalTransform(GeometryContext& gctx,
                                    const std::size_t portalIdx,
                                    Transform3&& portalLocToGlob) const = 0;

 private:
  std::vector<std::unique_ptr<detail::PortalPlacement>> m_portalPlacements{};
};
}  // namespace Acts
