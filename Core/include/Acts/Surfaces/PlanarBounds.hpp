// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <vector>

namespace Acts {

/// Forward declare rectangle bounds as boundary box
class RectangleBounds;

/// @class PlanarBounds
///
/// common base class for all bounds that are in a local x/y cartesian frame
///  - simply introduced to avoid wrong bound assigments to surfaces
///
class PlanarBounds : public SurfaceBounds {
 public:
  /// Return the vertices
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note that the extremas are given, which may slightly alter the
  /// number of segments returned
  ///
  /// @return vector for vertices in 2D
  virtual std::vector<Vector2D> vertices(unsigned int lseg = 1) const = 0;

  /// Mask a segment with the bounds shape, in case the segment is
  /// fully inside the bounds, no clipping is done. Both outside
  /// is currently not supported, it has to be checked upstream.
  ///
  /// @param start The start of the segment
  /// @param end The end of the segment
  ///
  /// @return a pair of the (potentially) clipped bounds
  virtual std::pair<Vector2D, Vector2D> mask(const Vector2D& start,
                                             const Vector2D& end) const;

  /// Bounding box parameters
  ///
  /// @return rectangle bounds for a bounding box
  virtual const RectangleBounds& boundingBox() const = 0;
};

inline std::pair<Vector2D, Vector2D> PlanarBounds::mask(
    const Vector2D& start, const Vector2D& end) const {
  bool startInside = inside(start, true);
  bool endInside = inside(end, true);
  if (not startInside or not endInside) {
    return detail::VerticesHelper::mask(vertices(1), start, startInside, end,
                                        endInside);
  }
  return {start, end};
}

}  // namespace Acts