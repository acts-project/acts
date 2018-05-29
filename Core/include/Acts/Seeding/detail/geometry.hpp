// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

namespace Acts {
namespace detail {

  /// Calculate the curvature of a circle from three points.
  ///
  /// @note Difference vectors must provide `.x()`, `.y()`
  ///
  /// Sign convention as defined in
  ///
  ///     P. Billoir and S. Qian, NIM A311 (1992), 139--150
  ///
  /// i.e. curvature is positive for a counter-clockwise rotation and
  /// negative for a clockwise rotation as defined by the order of the
  /// input points.
  template <typename D01, typename D12>
  inline double
  calcCircleCurvature(const D01& d01, const D12& d12)
  {
    double x01 = d01.x();
    double y01 = d01.y();
    double x12 = d12.x();
    double y12 = d12.y();
    double x02 = x01 + x12;
    double y02 = y01 + y12;
    // length of the triangle sides
    double a = std::hypot(x12, y12);
    double b = std::hypot(x02, y02);
    double c = std::hypot(x01, y01);
    // 2 * (signed) area of the triangle
    double k = (x02 * y01 - x01 * y02);
    // radius = product of side lengths / 4 times triangle area
    return (2 * k) / (a * b * c);
  }

  /// Estimate the distance of closest approach to the origin.
  ///
  /// @param p0     Initial circle position (must provide `.x()`, `.y()`)
  /// @param phi0   Initial azimuthal direction angle
  /// @param kappa  Signed circle curvature
  ///
  /// A first order approximation in the curvature of the analytical expression
  /// is used to support curvature -> 0.
  template <typename P0>
  inline double
  estimateCircleD0(const P0& p0, double phi0, double kappa)
  {
    double x0sin = p0.x() * std::sin(phi0);
    double y0cos = p0.y() * std::cos(phi0);
    // d0 = distance(origin, circle center) - radius
    // distance^2 / radius^2= (1 + 2*kappa*(y*cos(phi) - x*sin(phi))
    //                           + kappa^2*(x^2 + y^2))
    // approximate sqrt(distance^2 / radius^2) = sqrt(1 + x) = ...
    // TODO 2016-09-29 msmk: redid the calculations but distribution still weird
    return (y0cos - x0sin) + 0.5 * kappa * (p0.x() * p0.x() + p0.y() * p0.y())
        + 0.5 * kappa * (x0sin * x0sin + y0cos * y0cos - 2 * x0sin * y0cos);
  }

  /// Estimate the z position at vanishing transverse radius.
  template <typename P0>
  inline double
  estimateZHelixZ0(const P0& p0, double theta, double /*kappa*/)
  {
    return p0.z() - std::hypot(p0.x(), p0.y()) * std::tan(M_PI_2 - theta);
  }

  /// Calculate the intersection point of a line and a cylinder.
  ///
  /// @param p0      point on the line (must provide `.x()`, `.y()`, `.z()`)
  /// @param d       direction vector (must provide `.x()`, `.y()`, `.z()`)
  /// @param radius  radius of a circle centered around the origin
  ///
  /// Returns the intersection point with the closest propagation path
  /// along the line. If no intersection exists, the closest point on the line
  /// to the circle is returned. Intersection and closest approach are
  /// calculated only in the transverse plane but the final position
  /// is calculated in three dimensions.
  template <typename P0, typename D>
  inline Vector3D
  calcLineCircleIntersection(const P0& p0, const D& d, double radius)
  {
    // line equation x + d * u, with start point between the input points
    double x0 = p0.x();
    double y0 = p0.x();
    double z0 = p0.x();
    double dx = d.x();
    double dy = d.y();
    double dz = d.z();
    // intersection equation only in the transverse plane
    // (x + d * u)^2 - r^2 = 0
    // d^2 u^2 + 2 * x.d * u + x^2 - r^2 = alpha u^2 + beta u + gamma = 0
    double alpha = dx * dx + dy * dy;
    double beta  = 2 * (x0 * dx + y0 * dy);
    double gamma = x0 * x0 + y0 * y0 - radius * radius;
    double u     = 0;

    RealQuadraticEquation eq(alpha, beta, gamma);
    if (eq.solutions == 2) {
      u = (std::abs(eq.first) < std::abs(eq.second)) ? eq.first : eq.second;
    } else if (eq.solutions == 1) {
      u = eq.first;
    } else {
      // fall-back: closest approach to the circle, (x + d * u) . d = 0
      u = -(x0 * dx + y0 * dy) / (dx * dx + dy * dy);
    }
    return Vector3D(x0 + dx * u, y0 + dy * u, z0 + dz * u);
  }

}  // namespace detail
}  // namespace Acts