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

#include <ACTS/Surfaces/RealQuadraticEquation.hpp>
#include <ACTS/Utilities/Definitions.hpp>

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

    /**
     * @brief Calculate the intersection point of a line and a cylinder
     * @param l0 first point on the line (must provide `.x()`, `.y()`, `.z()`)
     * @param l1 second point on the line (must provide `.x()`, `.y()`, `.z()`)
     * @param radius radius of a circle centered around the origin
     *
     * Returns the intersection point with the closest propagation path
     * along the line. If no intersection exists, the closest point on the line
     * to the circle is returned. Intersection and closest approach are
     * calculated only in the transverse plane but the final position
     * is calculated in three dimensions.
     */
    template <typename L0, typename L1>
    Acts::Vector3D
    calcLineCircleIntersection(const L0& l0, const L1& l1, double radius)
    {
      // line equation x + d * u, with start point between the input points
      auto x0 = (l0.x() + l1.x()) / 2;
      auto y0 = (l0.y() + l1.y()) / 2;
      auto z0 = (l0.z() + l1.z()) / 2;
      auto dx = l1.x() - x0;
      auto dy = l1.y() - y0;
      auto dz = l1.z() - z0;
      // intersection equation only in the transverse plane
      // (x + d * u)^2 - r^2 = 0
      // d^2 u^2 + 2 * x.d * u + x^2 - r^2 = alpha u^2 + beta u + gamma = 0
      double alpha = dx * dx + dy * dy;
      double beta  = 2 * (x0 * dx + y0 * dy);
      double gamma = x0 * x0 + y0 * y0 - radius * radius;
      double u     = 0;

      RealQuadraticEquation eq(alpha, beta, gamma);
      if (eq.solutions = two) {
        u = (std::abs(eq.first) < std::abs(eq.second)) ? eq.first : eq.second;
      } else if (eq.solutions = one) {
        u = eq.first;
      } else {
        // fall-back: closest approach to the circle, (x + d * u) . d = 0
        u = -(x0 * dx + y0 * dy) / (dx * dx + dy * dy);
      }
      return Acts::Vector3D(x0 + dx * u, y0 + dy * u, z0 + dz * u);
    }

  }  // namespace detail
}  // namespace Seeding
}  // namespace Acts

#endif  // ACTS_SEEDING_GEOMETRY_HPP
