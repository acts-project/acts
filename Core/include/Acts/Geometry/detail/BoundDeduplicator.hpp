// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometryVisitor.hpp"
#include "Acts/Utilities/BoundFactory.hpp"

namespace Acts::detail {
/// @brief Tracking geometry visitor that deduplicates the bounds of Surfaces &
///        TrackingVolumes. E.g., if two PlaneSurfaces have each a
///        RectangularBounds object with (10.cm, 5.cm), the visitor ensures that
///        there is only one instance of with these parameters and that both
///        surfaces hold this pointer. The same logic applies to the sharing
///        of the TrackingVolumeBounds
class BoundDeduplicator : public TrackingGeometryMutableVisitor {
 public:
  /// @brief Visit and potentially modify a tracking volume
  /// @param volume The tracking volume being visited
  /// @note Called for each volume in the geometry hierarchy during traversal
  void visitVolume(TrackingVolume& volume) final;

  /// @brief Visit and potentially modify a portal
  /// @param portal The portal being visited
  /// @note Called for each portal encountered during geometry traversal
  void visitPortal(Portal& portal) final;

  /// @brief Visit and potentially modify a surface
  /// @param surface The surface being visited
  /// @note Called for each surface encountered during geometry traversal
  void visitSurface(Surface& surface) final;
  /// @brief Visit and potentially modify a boundary surface
  /// @param boundary The boundary surface being visited
  /// @note Called for each boundary surface encountered during geometry traversal
  void visitBoundarySurface(BoundarySurfaceT<TrackingVolume>& boundary) final;

 private:
  SurfaceBoundFactory m_surfFactory;
  VolumeBoundFactory m_volFactory;
};
}  // namespace Acts::detail
