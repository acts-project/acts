// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/detail/CurvatureCorrection.hpp"

namespace Acts::detail::CylinderHelper {

// Inputs are in the cylinder's local frame. Preserve arrival at the actual
// starting position; only predicted crossings receive a correction.
inline std::optional<std::pair<Vector3, double>> correctIntersection(
    double s, const Vector3& seed, const Vector3& direction,
    const Vector3& curvature, double curvatureNorm, double radius,
    double tolerance) {
  if (std::abs(s) < std::abs(tolerance)) {
    return std::nullopt;
  }
  // n = (seed.x, seed.y, 0)/R. The common 1/R cancels in ds; using
  // the actual seed avoids cancellation from sharing expanded dot products
  // between distant ray roots. The incidence guard is scaled by R too.
  return correctCurvatureIntersection(
      s, seed, direction, curvature, seed.head<2>().dot(direction.head<2>()),
      seed.head<2>().dot(curvature.head<2>()), curvatureNorm, radius);
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

}  // namespace Acts::detail::CylinderHelper
