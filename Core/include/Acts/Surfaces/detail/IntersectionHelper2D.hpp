// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <utility>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {
namespace detail {

struct IntersectionHelper2D {
  /// Intersect two segments
  ///
  /// @param s0 The Start of the segement
  /// @param s1 The end of the segement
  /// @param origin The Start of intersection line
  /// @param direction The Direction of intersection line
  ///
  /// @return the intersection point with status
  static Intersection2D intersectSegment(const Vector2D& s0, const Vector2D& s1,
                                         const Vector2D& origin,
                                         const Vector2D& dir);

  /// Intersect ellipses
  ///
  /// @param Rx The radius in x
  /// @param Ry The radius in y
  /// @param origin The Start of intersection line
  /// @param direction The Direction of intersection line
  ///
  /// @return the intersection points
  static std::pair<Intersection2D, Intersection2D> intersectEllipse(
      double Rx, double Ry, const Vector2D& origin, const Vector2D& dir);

  /// Intersect the circle
  ///
  /// @param R The radius
  /// @param origin The Start of intersection line
  /// @param direction The Direction of intersection line
  ///
  /// @return the intersection points
  static inline std::pair<Intersection2D, Intersection2D> intersectCircle(
      double R, const Vector2D& origin, const Vector2D& dir) {
    return intersectEllipse(R, R, origin, dir);
  }

};  // struct IntersectionHelper2D

}  // namespace detail
}  // namespace Acts
