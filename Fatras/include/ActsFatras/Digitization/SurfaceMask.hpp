// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/detail/IntersectionHelper2D.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <vector>

namespace Acts {
class Surface;
class AnnulusBounds;
class RadialBounds;
class CylinderBounds;
}  // namespace Acts

namespace ActsFatras {

/// A brief struct that allows to apply a surface bound mask.
///
/// A single implementation handles all supported surface types, dispatched
/// internally on `surface.type()`:
///   - Plane / Disc / DiscTrapezoid : polygon, radial or annulus masking in
///     the Cartesian / polar local frame.
///   - Cylinder : axis-aligned rectangle clipping in the unrolled (rPhi, z)
///     readout frame.
struct SurfaceMask {
  /// Shorthand for a 2-d segment;
  using Segment2D = std::array<Acts::Vector2, 2>;

  /// Apply the mask on the segment
  /// - If the segment is fully inside the surface, return unchanged
  /// - Otherwise mask/clip the segment to fit into the bounds
  ///
  /// @note PlaneSurface/DiscSurface/CylinderSurface are supported
  ///
  /// @note If both end points of the segment are inside, the segment
  /// is not clipped/masked, even if it would cross a surface boundary.
  /// Examples for those would be non-covex polygons or segments on a
  /// radial bound, where the radial boundary is crossed. Such segments
  /// do not occur in Digitization, as the hit has to be inside the
  /// surface bounds to start with.
  ///
  /// @param surface The surface in question
  /// @param segment The track segment (on surface)
  ///
  /// @return a result wrapping a segment
  Acts::Result<Segment2D> apply(const Acts::Surface& surface,
                                const Segment2D& segment) const;

  /// Apply the mask of a polygon
  ///
  /// @param vertices The vertices of the polygon
  /// @param segment The track segment (on surface)
  /// @param firstInside The indicator if the first is inside
  ///
  /// @return a result wrapping a segment
  Acts::Result<Segment2D> polygonMask(
      const std::vector<Acts::Vector2>& vertices, const Segment2D& segment,
      bool firstInside) const;

  /// Apply the mask of a Radial disk
  ///
  /// @param rBounds The radial disc for the masking
  /// @param segment The track segment (on surface)
  /// @param polarSegment The track segment (on surface, in polar)
  /// @param firstInside The indicator if the first is inside
  ///
  /// @return a result wrapping a segment
  Acts::Result<Segment2D> radialMask(const Acts::RadialBounds& rBounds,
                                     const Segment2D& segment,
                                     const Segment2D& polarSegment,
                                     bool firstInside) const;

  /// Apply the mask of an annulus disk
  ///
  /// @param aBounds The annulus disc for the masking
  /// @param segment The track segment (on surface)
  /// @param firstInside The indicator if the first is inside
  ///
  /// @return a result wrapping a segment
  Acts::Result<Segment2D> annulusMask(const Acts::AnnulusBounds& aBounds,
                                      const Segment2D& segment,
                                      bool firstInside) const;

  /// Apply the mask of a cylinder (axis-aligned rectangle in (rPhi, z))
  ///
  /// @param cBounds The cylinder bounds for the masking
  /// @param segment The track segment in the unrolled (rPhi, z) frame
  ///
  /// @return a result wrapping a segment
  Acts::Result<Segment2D> cylinderMask(const Acts::CylinderBounds& cBounds,
                                       const Segment2D& segment) const;

  /// 2D intersection helper for geometric calculations
  Acts::detail::IntersectionHelper2D intersector{};
};

}  // namespace ActsFatras
