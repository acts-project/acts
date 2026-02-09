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

/// @brief Interface class to define the transform cache backend for alignable volumes.
///        Alignable volumes can be dynamically moved by the client code using
///        the information wrapped into the GeometryContext similar to the
///        alignable surfaces. The challenge is that the boundary surfaces of
///        the volume need to move accordingly and that ACTS does not know
///        anything about the caching of the alignable geometry objects a
///        priori.
///
///        A client-based implementation of the VolumePlacementBase can be
///        parsed to the constructor of the Volume instead of the fixed
///        transform. The Volume is then querying its global position from the
///        VolumePlacements. The associated bonudary surfaces are also
///        requesting their position in global space. The Volume's transform and
///        the bounds can the no longer be overwritten at a later stage.
///
///        An implementation of the VolumePlacementBase needs to satisfy the
///        following interface.
///
///         1) Transforms switching from the volume's frame into the global
///            experiment's frame and vice versa:
///
///           const Transform3& localToGlobalTransform(const GeometryContext&
///           gctx) const;
///
///           const Transform3& localToGlobalTransform(const GeometryContext&
///           gctx) const;
///
///
///        2) At the end of the tracking geometry construction, the portals
///           that are associated to the volume aligned by the
///           VolumePlacementBase are connected to this particular instance.
///           Via the
///
///              void expandTransformCache(const GeometryContext& gctx,
///                                        const std::size_t nPortals);
///
///           method, ACTS requests the client to expand the cache backend
///           by nPortals where the associated portal transforms are cached
///           and updated with the change of the alignment constants. At
///           creation, the cache needs to be filled with the initial transforms
///           as the framework may check that the portal did not move during the
///           procedure
///
///             const Transform3& portalLocalToGlobal(const GeometryContext&
///             gctx,
///                                                   const std::size_t
///                                                   portalIdx) const;
///
///           Is called every time when the portal is asking for its position in
///           space. The portalIdx is the unique index of the portal and assists
///           the client to return the appropriate transform
///
///       3)  Every time when the alignment of the volume is updated, the client
///           also needs to cache the transforms of the associated surfaces.
///           After, the central volume has moved, a loop similar to the one
///           below needs to be implemented:
///
///             for (std::size_t p = 0 ; p < nPortalPlacements(); ++p) {
///                 context.cacheInContext(context.getContext(), p);
///             }
///             context.volGlobToLocal = context.volLocToGlobal.inverse();
///

class VolumePlacementBase {
 public:
  /// @brief Default constructor
  VolumePlacementBase() noexcept;

  /// @brief Virtual default destructor
  virtual ~VolumePlacementBase();

  /// @brief Move constructor
  VolumePlacementBase(VolumePlacementBase&& other) noexcept;

  /// @brief Move assignment operator
  VolumePlacementBase& operator=(VolumePlacementBase&& other) noexcept;

  /// @brief Delete copy constructor
  VolumePlacementBase(const VolumePlacementBase& other) = delete;

  /// @brief Delete copy assignment operator
  VolumePlacementBase& operator=(const VolumePlacementBase& other) = delete;

  /// @brief Abrivation of the portal surface vector
  using PortalVec_t = std::vector<std::shared_ptr<RegularSurface>>;

  /// @brief Receives the vector of oriented portal surfaces produced by the
  ///        VolumeBounds and makes them to float with the alignment provided
  ///        by the volume. It then the vector of updated oriented surfaces
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param portalsToAlign: List of portals to align
  void makePortalsAlignable(const GeometryContext& gctx,
                            const PortalVec_t& portalsToAlign);

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

  /// @brief Returns the const pointer to the `SurfacePlacementBase` object
  ///        aligning the i-th portal (May be nullptr if index exceeds the
  ///        number of portals)
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  const detail::PortalPlacement* portalPlacement(
      const std::size_t portalIdx) const;

 protected:
  /// @brief This method is called when the portal surfaces are registered with
  ///        the VolumePlacmentBase class. The client receives the information
  ///        about how many portals exist and he is requested to adapt the size
  ///        of the transform cache backend accordingly. The method is called
  ///        after the portal placements have been created.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param nPortals: The number of portals that are registered with this volume
  ///                  placement instance
  virtual void expandTransformCache(const GeometryContext& gctx,
                                    const std::size_t nPortals) = 0;

  /// @brief Returns the transform from the portal's frame into the
  ///        experiment's global frame taking the alignment corrections
  ///        of the associated volume into account.
  ///        @note: The call of this function is only allowed after the
  ///               volume itself is moved. A swapped call order probably
  ///               leads to unaligned portals
  /// @param gctx: The geometry context carrying the current volume alignment
  /// @param portalIdx: Index of the portal to align
  Transform3 alignPortal(const GeometryContext& gctx,
                         const std::size_t portalIdx) const;

 private:
  std::vector<std::unique_ptr<detail::PortalPlacement>> m_portalPlacements{};
};
}  // namespace Acts
