// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/detail/CurvatureCorrection.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

namespace Acts::detail::CylinderHelper {

/// Solve |position + s*direction|^2 = radius^2 in the transverse plane.
/// Inputs may be local xy vectors or global vectors crossed with the cylinder
/// axis. Keeping that choice with the caller avoids redundant transformations.
template <typename PositionDerived, typename DirectionDerived>
inline RealQuadraticEquation intersectionSolver(
    const Eigen::MatrixBase<PositionDerived>& position,
    const Eigen::MatrixBase<DirectionDerived>& direction, double radius) {
  return RealQuadraticEquation(direction.squaredNorm(),
                               2. * direction.dot(position),
                               position.squaredNorm() - radius * radius);
}

/// Inputs are in the cylinder's local frame. Preserve arrival at the actual
/// starting position; only predicted crossings receive a correction.
inline std::optional<CurvatureCorrection> correctIntersection(
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

/// Keep projection policy internal: subsequent numerical optimizations can
/// replace this implementation without changing Surface or navigation APIs.
/// Corrected cylinder points lie on the seed's tangent plane. Radially project
/// them back onto the cylinder, leaving z and the estimated path unchanged.
/// A bound-coordinate Jacobian could instead map the tangential displacement
/// u to an arc-length increment u. Radial projection gives R*atan(u/R), so the
/// two agree to first order but differ beyond it. Keep the Cartesian correction
/// and radial projection together to preserve the current approximation.
/// This does not make the trajectory approximation second-order accurate.
struct CylinderProjectionScales {
  double first = 0.;
  double second = 0.;
};

/// Exact paired radial scales; zero marks an invalid projection.
inline CylinderProjectionScales cylinderProjectionScales(
    const std::optional<CurvatureCorrection>& first,
    const std::optional<CurvatureCorrection>& second, double radius) {
  // Calculate both independent lengths before dividing. This allows the
  // compiler to overlap the two exact projections, without SIMD or an
  // approximate reciprocal square root. Zero marks a failed projection.
  const double r1 = first ? first->position.head<2>().norm() : 0.;
  const double r2 = second ? second->position.head<2>().norm() : 0.;
  return {r1 > 0. && std::isfinite(r1) ? radius / r1 : 0.,
          r2 > 0. && std::isfinite(r2) ? radius / r2 : 0.};
}

}  // namespace Acts::detail::CylinderHelper
