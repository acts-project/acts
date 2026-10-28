// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <cmath>
#include <optional>
#include <utility>

namespace Acts::detail {

// Local trajectory model: x(u) = x0 + u*t + u^2*k/2, with unit t and
// geometric curvature k = d^2 x / ds^2. At the ray seed y = x0 + s*t,
// write the corrected path as s + ds. To first order in k,
//
//   x(s + ds) = y + (s^2/2)*k + ds*t + O(k^2),
//   n.(x(s + ds) - y) = 0,
//   ds = -(s^2/2)*(n.k)/(n.t).
//
// This drops s*ds*k and ds^2*k/2. It is a tangent-plane perturbation, not
// an exact helix intersection. The incidence, bend and relative path-shift
// guards below limit the perturbation; they do NOT bound its physical error.
// A scaled normal is allowed: normalScale is its length. For cylinders the
// radial seed vector can be used directly, cancelling the common 1/R factor
// in the ratio. All vectors must use the same coordinate frame.
inline std::optional<std::pair<Vector3, double>> correctCurvatureIntersection(
    double s, const Vector3& seed, const Vector3& direction,
    const Vector3& curvature, double incidence, double normalCurvature,
    double curvatureNorm, double normalScale = 1.) {
  if (std::abs(incidence) < 1e-6 * normalScale ||
      std::abs(s) * curvatureNorm > 0.2) {
    return std::nullopt;
  }
  const double halfPathSquared = 0.5 * s * s;
  const double ds = -halfPathSquared * normalCurvature / incidence;
  if (!std::isfinite(ds) || std::abs(ds) > 0.2 * std::abs(s)) {
    return std::nullopt;
  }
  const Vector3 point = seed + halfPathSquared * curvature + ds * direction;
  if (!point.allFinite()) {
    return std::nullopt;
  }
  return std::pair<Vector3, double>{point, s + ds};
}

inline Intersection3D correctPlanarCurvatureIntersection(
    const Intersection3D& seed, const Vector3& direction,
    const Vector3& curvature, const Vector3& normal) {
  // Arrival is determined from the actual starting position, never predicted.
  if (!seed.isValid() || seed.status() == IntersectionStatus::onSurface) {
    return seed;
  }
  if (auto corrected = correctCurvatureIntersection(
          seed.pathLength(), seed.position(), direction, curvature,
          normal.dot(direction), normal.dot(curvature), curvature.norm())) {
    return {corrected->first, corrected->second, IntersectionStatus::reachable};
  }
  return seed;
}

// Keep projection policy internal: subsequent numerical optimizations can
// replace this implementation without changing Surface or navigation APIs.
// Corrected cylinder points lie on the seed's tangent plane. Radially project
// them back onto the cylinder, leaving z and the estimated path unchanged.
// This does not make the trajectory approximation second-order accurate.
struct CylinderProjectionScales {
  double first = 0.;
  double second = 0.;
};

inline CylinderProjectionScales cylinderProjectionScales(
    const std::optional<std::pair<Vector3, double>>& first,
    const std::optional<std::pair<Vector3, double>>& second, double radius) {
  // Calculate both independent lengths before dividing. This allows the
  // compiler to overlap the two exact projections, without SIMD or an
  // approximate reciprocal square root. Zero marks a failed projection.
  const double r1 = first ? first->first.head<2>().norm() : 0.;
  const double r2 = second ? second->first.head<2>().norm() : 0.;
  return {r1 > 0. && std::isfinite(r1) ? radius / r1 : 0.,
          r2 > 0. && std::isfinite(r2) ? radius / r2 : 0.};
}

}  // namespace Acts::detail
