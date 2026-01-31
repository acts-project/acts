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
struct OrientedSurface;

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
  VolumePlacementBase() noexcept;

  /// @brief Move constructor
  VolumePlacementBase(VolumePlacementBase&& other) noexcept;

  /// @brief Move assignment operator
  VolumePlacementBase& operator=(VolumePlacementBase&& other) noexcept;

  /// @brief Default destructor
  virtual ~VolumePlacementBase() = default;

  /// @brief Abrivation of the portal surface vector
  using PortalVec_t = std::vector<OrientedSurface>;

  /// @brief Receives the vector of oriented portal surfaces produced by the
  ///        VolumeBounds and makes them to float with the alignment provided
  ///        by the volume. It then the vector of updated oriented surfaces
  /// @param portalsToAlign: List of portals to align
  PortalVec_t makePortalsAlignable(PortalVec_t&& portalsToAlign);

  /// @brief Returns the number of portal placement objects
  std::size_t nPortalPlacements() const;

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

 protected:
  /// @brief Returns the const pointer to the `SurfacePlacementBase` object
  ///        aligning the i-th portal (May be nullptr if index exceeds the
  ///        number of portals)
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  const detail::PortalPlacement* portalPlacement(
      const std::size_t portalIdx) const;

  /// @brief This method is called for the first time when the portal surfaces
  ///        are registered with the VolumePlacmentBase class. The client
  ///        receives the information about how many portals exist and it's
  ///        requested to adapt the size of the transform cache backend
  ///        accordingly. The method is called after the portal placements have
  ///        been created.
  /// @param nPortals: The number of portals that are registered with this volume
  ///                  placement instance
  virtual void expandTransformCache(const std::size_t nPortals) = 0;

  /// @brief Returns the
  Transform3 alignPortal(const GeometryContext& gctx,
                         const std::size_t portalIdx) const;

 private:
  std::vector<std::unique_ptr<detail::PortalPlacement>> m_portalPlacements{};
};
}  // namespace Acts
