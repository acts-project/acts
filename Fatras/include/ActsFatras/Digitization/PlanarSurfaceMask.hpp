// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>

namespace Acts {
class Surface;
}

namespace ActsFatras {

/// A brief struct that allows to apply a surface bound mask.
struct PlanarSurfaceMask {
  /// Shorthand for a 2-d segment;
  using Segment2D = std::array<Acts::Vector2D, 2>;

  /// Apply the mask on the segment
  /// - If the semgent is full inside the surface, return unchanged
  /// - Otherwise mask/clip the segment to fit into the bounds
  ///
  /// @note Only PlaneSurface/DiscSurface are supported
  ///
  /// @param geoCtx The geometry context
  /// @param surface The surface in question
  /// @param segment The track segment (on surface)
  ///
  /// @return a result wrapping a segment
  Acts::Result<Segment2D> apply(const Acts::Surface& surface,
                                const Segment2D& segment) const;

  /// Apply the mask of a polygon
  ///
  /// @param outise The outside point of the segment
  /// @param inside The inside point of the segment
  /// @param vertices The vertices of the polygon
  ///
  /// @return a result wrapping the new new outside position
  Acts::Result<Acts::Vector2D> polygonMask(
      const Acts::Vector2D& outside, const Acts::Vector2D& inside,
      const std::vector<Acts::Vector2D>& vertices) const;
};

struct PolarSegment {
  Acts::Vector2D clipped;
  Acts::Vector2D start;
  Acts::Vector2D end;
  bool inside = false;

  PolarSegment(Acts::Vector2D clipped_, const Acts::Vector2D& start_,
               const Acts::Vector2D& end_, bool inside_)
      : clipped(clipped_), start(start_), end(end_), inside(inside_) {}
};

}  // namespace ActsFatras
