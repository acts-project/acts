// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"

namespace Acts {

class TrackingVolume;
class GridPortalLink;

/// Trivial portal link links to a single target volume on every point on a
/// surface
class TrivialPortalLink final : public PortalLinkBase {
 public:
  /// Construct a trivial portal link from a surface and a volume
  /// @param surface is the surface
  /// @param volume is the target
  TrivialPortalLink(std::shared_ptr<RegularSurface> surface,
                    TrackingVolume& volume)
      : PortalLinkBase(std::move(surface)), m_volume{&volume} {}

  /// Make a 1D grid portal link from this trivial portal link
  /// The grid size is automatically determined from the surface bounds.
  /// @param direction The binning direction
  /// @return A grid
  std::unique_ptr<GridPortalLink> makeGrid(BinningValue direction) const;

  /// Print the portal link to a stream
  /// @param os output stream
  void toStream(std::ostream& os) const override;

  /// Resolve the volume for a 2D position
  /// @note Always returns the single target volume
  /// @param gctx is the geometry context
  /// @param position is the 2D position
  /// @param tolerance is the tolerance
  /// @return The target volume (can be null)
  Result<const TrackingVolume*> resolveVolume(
      const GeometryContext& gctx, const Vector2& position,
      double tolerance = s_onSurfaceTolerance) const override;

  /// Resolve the volume for a 3D position
  /// @note Always returns the single target volume
  /// @param gctx is the geometry context
  /// @param position is the 2D position
  /// @param tolerance is the tolerance
  /// @return The target volume (can be null)
  /// @note The position is assumed to be on the associated surface.
  Result<const TrackingVolume*> resolveVolume(
      const GeometryContext& gctx, const Vector3& position,
      double tolerance = s_onSurfaceTolerance) const override;

  /// Get the single volume that this trivial portal link is associated with
  /// @return The target volume
  const TrackingVolume& volume() const;

 private:
  TrackingVolume* m_volume;
};

}  // namespace Acts
