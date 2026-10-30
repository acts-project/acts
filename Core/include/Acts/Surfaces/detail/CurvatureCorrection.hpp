// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>
#include <optional>

namespace Acts::detail {

/// Predicted position and signed path before shape projection and bounds
/// checks.
struct CurvatureCorrection {
  /// Position before any cylinder projection.
  Vector3 position;
  /// Signed corrected path from the original starting point.
  double pathLength;
};

/// Reject near-tangent seeds: below this cosine, the normal equation can
/// amplify displacement errors by more than 10^6. This is a conditioning guard.
static constexpr double s_minCurvatureIncidence = 1e-6;
/// Limit the estimated tangent change to 0.2 (~11 degrees for a circle).
/// This keeps the frozen-curvature expansion local; it is not an error bound.
static constexpr double s_maxCurvatureBend = 0.2;
/// Limit the path adjustment to 20% of the ray path, keeping the prediction
/// near its seed branch. This is a heuristic, not a calibrated accuracy cut.
static constexpr double s_maxRelativeCurvaturePathShift = 0.2;

/// Compile-time selection for all curvature-aware surface intersections (and
/// hence navigation). Changing this requires rebuilding the library.
static constexpr bool s_useQuadraticCurvatureCorrection = false;

/// Local trajectory model: x(u) = x0 + u*t + u^2*k/2, with unit t and
/// geometric curvature k = d^2 x / ds^2. The vector t is the starting track
/// direction, not a direction along the surface. The total shift
/// (s^2/2)*k + ds*t lies in the surface tangent plane. At y = x0 + s*t,
/// write the corrected path as s + ds. To first order in k,
///
///   x(s + ds) = y + (s^2/2)*k + ds*t + O(k^2),
///   n.(x(s + ds) - y) = 0,
///   ds = -(s^2/2)*(n.k)/(n.t).
///
/// This drops s*ds*k and ds^2*k/2. It is a tangent-plane perturbation, not
/// an exact helix intersection. The incidence, bend and relative path-shift
/// guards below limit the perturbation; they do NOT bound its physical error.
/// A scaled normal is allowed: normalScale is its length. For cylinders the
/// radial seed vector can be used directly, cancelling the common 1/R factor
/// in the ratio. All vectors must use the same coordinate frame.
/// @tparam Quadratic Retain all local-parabola terms in the tangent-plane solve
/// @param rayPathLength Signed ray path from the starting position to the seed
/// @param rayPosition Ray intersection position, in the same frame as all vectors
/// @param direction Unit track direction at the starting position
/// @param curvature Starting geometric curvature d^2 x / ds^2 (inverse length)
/// @param incidence n.dot(direction), using the possibly scaled surface normal n
/// @param normalCurvature n.dot(curvature), using the same normal
/// @param curvatureNorm Magnitude of curvature, reusable across both ray roots
/// @param normalScale Magnitude of n (one for a unit normal)
/// @return Prediction, or nullopt when the caller should retain the ray seed
template <bool Quadratic = s_useQuadraticCurvatureCorrection>
inline std::optional<CurvatureCorrection> correctCurvatureIntersection(
    double rayPathLength, const Vector3& rayPosition, const Vector3& direction,
    const Vector3& curvature, double incidence, double normalCurvature,
    double curvatureNorm, double normalScale = 1.) {
  if (std::abs(incidence) < s_minCurvatureIncidence * normalScale ||
      std::abs(rayPathLength) * curvatureNorm > s_maxCurvatureBend) {
    return std::nullopt;
  }
  const double halfPathSquared = 0.5 * rayPathLength * rayPathLength;
  double ds = 0.;
  double displacementScale = halfPathSquared;
  if constexpr (Quadratic) {
    // With s = rayPathLength, retain s*ds*k and ds^2*k/2 in the parabola:
    //   a*ds^2 + b*ds + c = 0,
    //   a = n.k/2, b = n.t + s*n.k, c = (s^2/2)*n.k.
    // The rationalized root is continuous with ds = 0 as k -> 0.
    // A reversed slope or negative discriminant invalidates this branch.
    // This solves the seed tangent plane, not the full cylinder or helix.
    const double a = 0.5 * normalCurvature;
    const double b = incidence + rayPathLength * normalCurvature;
    const double c = halfPathSquared * normalCurvature;
    const double discriminant = b * b - 4. * a * c;
    if (discriminant < 0. || b * incidence <= 0.) {
      return std::nullopt;
    }
    ds = -2. * c / (b + std::copysign(std::sqrt(discriminant), b));
    displacementScale = 0.5 * (rayPathLength + ds) * (rayPathLength + ds);
  } else {
    ds = -halfPathSquared * normalCurvature / incidence;
  }
  if (!std::isfinite(ds) || std::abs(ds) > s_maxRelativeCurvaturePathShift *
                                               std::abs(rayPathLength)) {
    return std::nullopt;
  }
  const Vector3 point =
      rayPosition + displacementScale * curvature + ds * direction;
  if (!point.allFinite()) {
    return std::nullopt;
  }
  return CurvatureCorrection{point, rayPathLength + ds};
}

}  // namespace Acts::detail
