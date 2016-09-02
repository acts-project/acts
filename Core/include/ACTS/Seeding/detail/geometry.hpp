// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SEEDING_GEOMETRY_HPP
#define ACTS_SEEDING_GEOMETRY_HPP

#include <cmath>

namespace Acts {
namespace Seeding {
  namespace detail {

    /**
     * @brief Calculate the curvature of a circle from three points.
     */
    template <typename P0, typename P1, typename P2>
    double
    calcCircleCurvature(const P0& p0, const P1& p1, const P2& p2)
    {
      double x01 = p1.x() - p0.x();
      double y01 = p1.y() - p0.y();
      double x02 = p2.x() - p0.x();
      double y02 = p2.y() - p0.y();
      // length of the triangle sides
      double a = std::hypot(p2.x() - p1.x(), p2.y() - p1.y());
      double b = std::hypot(x02, y02);
      double c = std::hypot(x01, y01);
      // 2 * (signed) area of the triangle
      double k = (x01 * y02 - x02 * y01);
      // radius = product of side lengths / 4 times triangle area
      // TODO check sign convention, clockwise-/counterclockwise
      return (2 * k) / (a * b * c);
    }

  }  // namespace detail
}  // namespace Seeding
}  // namespace Acts

#endif  // ACTS_SEEDING_GEOMETRY_HPP
