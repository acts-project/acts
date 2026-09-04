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

/// Interface class to define the transform cache backend for alignable volumes.
/// Alignable volumes can be dynamically moved by the client code using
/// the information wrapped into the GeometryContext similar to the
/// alignable surfaces. The challenge is that the boundary surfaces of
/// the volume need to move accordingly and that ACTS does not know
/// anything about the caching of the alignable geometry objects a
/// priori.
///
/// A client-based implementation of the VolumePlacementBase can be
/// passed to the constructor of the Volume instead of the fixed
/// transform. The Volume is then querying its global position from the
/// VolumePlacements. The associated bonudary surfaces are also
/// requesting their position in global space. The Volume's transform and
/// the bounds can the no longer be overwritten at a later stage.
///
/// An implementation of the @ref VolumePlacementBase needs to satisfy the
/// following interface.
///
///  1) Transforms switching from the volume's frame into the global
///     experiment's frame and vice versa:
///
///      const Transform3& localToGlobalTransform(const GeometryContext& gctx)
///      const;
///
///      const Transform3& localToGlobalTransform(const GeometryContext& gctx)
///      const;
///
///
///  2) At the end of the tracking geometry construction, the portals
///     that are associated to the volume aligned by the
///     VolumePlacementBase are connected to this particular instance.
///     The user may override the @ref makePortalsAlignable method to
///     instantiate a customized transform cache backend. He needs to ensure
///     that the base definition of the method is called from his implementation
///     otherwise the placements aligning the portals are note created.
///
///     Every time when the portal is asked for its position in space,
///     it's forwarding the request to the VolumePlacementBase by calling the
///
///         const Transform3& portalLocalToGlobal(const GeometryContext& gctx,
///                                               const std::size_t portalIdx)
///                                               const;
///
///     method. The portalIdx is the unique index of the portal and assists
///     the client to return the appropriate transform
///
///  3)  Every time when the alignment of the volume is updated, the client
///      also needs to cache the transforms of the associated surfaces.
///      After, the central volume has moved, a loop similar to the one
///      below needs to be implemented:
///
///        for (std::size_t p = 0 ; p < nPortalPlacements(); ++p) {
///            context.cachePortal(alignPortal(Acts::Geometrycontext{context},
///            p), p);
///        }
///        context.volGlobToLocal = context.volLocToGlobal.inverse();
///
///      The @ref alignPortal is a wrapper method attaching the portal -> volume transform
///      to the aligned local -> global transform of the volume and returning
///      the result via copy.
class VolumePlacementBase {
 public:
  /// Default constructor
  VolumePlacementBase() noexcept;

  /// Virtual default destructor
  virtual ~VolumePlacementBase();

  /// Delete the move constructor
  VolumePlacementBase(VolumePlacementBase&&) = delete;
  /// Delete the copy constructor
  VolumePlacementBase(const VolumePlacementBase&) = delete;

  /// Delete move assignment
  VolumePlacementBase& operator=(VolumePlacementBase&&) = delete;
  /// Delete copy assignment
  VolumePlacementBase& operator=(const VolumePlacementBase&) = delete;

  /// Receives the vector of oriented portal surfaces produced by the
  /// VolumeBounds and makes them to float with the alignment provided
  /// by the volume. It then the vector of updated oriented surfaces
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param portalsToAlign: List of portals to align
  virtual void makePortalsAlignable(
      const GeometryContext& gctx,
      const std::vector<std::shared_ptr<RegularSurface>>& portalsToAlign);

  /// Number of registered SurfacePlacement objects aligning the
  /// associated portals with the volume
  /// @returns The number of registered portal placements
  std::size_t nPortalPlacements() const;

  /// Returns the transformation from the local volume coordinates to
  ///        the experiment's global coordinate system
  /// @param gctx The current geometry context object, e.g. alignment
  /// @returns Reference to the local -> global transform
  virtual const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const = 0;

  /// Returns the transformation from the experiment's global frame to the
  /// local volume coordinate system
  /// @param gctx The current geometry context object, e.g. alignment
  /// @returns Reference to the global -> local transform
  virtual const Transform3& globalToLocalTransform(
      const GeometryContext& gctx) const = 0;

  /// Returns the transform from the portal's frame to the experiment's
  /// global frame for the portal surface associated with the volume
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  /// @returns Reference to the local -> global transform of the i-th portal
  virtual const Transform3& portalLocalToGlobal(
      const GeometryContext& gctx, const std::size_t portalIdx) const = 0;

  /// Pointer to the `SurfacePlacementBase` object aligning the i-th portal
  /// (May be nullptr if index exceeds the number of portals)
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  /// @returns Pointer to the i-th portal placement
  const detail::PortalPlacement* portalPlacement(
      const std::size_t portalIdx) const;
  /// Pointer to the `SurfacePlacementBase` object aligning the i-th portal
  /// (May be nullptr if index exceeds the number of portals)
  /// @param portalIdx: Internal index of the portal surface [0 - number of portals)
  /// @returns Pointer to the i-th portal placement
  detail::PortalPlacement* portalPlacement(const std::size_t portalIdx);

 protected:
  /// Constructs the transform from the portal's frame into the
  /// experiment's global frame taking the alignment corrections
  /// of the associated volume into account.
  /// @note: The call of this function is only allowed after the
  ///        volume itself is moved. A swapped call order probably
  ///        leads to unaligned portals
  /// @param gctx: The geometry context carrying the current volume alignment
  /// @param portalIdx: Index of the portal to align
  /// @returns The aligned localToGlobalTransform of the i-the portal
  Transform3 alignPortal(const GeometryContext& gctx,
                         const std::size_t portalIdx) const;

 private:
  /// Resource allocation of the SurfacePlacements to align the portals
  std::vector<std::unique_ptr<detail::PortalPlacement>> m_portalPlacements{};
};
}  // namespace Acts
