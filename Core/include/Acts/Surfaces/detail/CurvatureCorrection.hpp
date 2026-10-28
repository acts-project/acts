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

// Compile-time selection for all curvature-aware surface intersections (and
// hence navigation). Changing this requires rebuilding the library.
static constexpr bool s_useQuadraticCurvatureCorrection = false;

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
template <bool Quadratic = s_useQuadraticCurvatureCorrection>
inline std::optional<std::pair<Vector3, double>> correctCurvatureIntersection(
    double s, const Vector3& seed, const Vector3& direction,
    const Vector3& curvature, double incidence, double normalCurvature,
    double curvatureNorm, double normalScale = 1.) {
  if (std::abs(incidence) < 1e-6 * normalScale ||
      std::abs(s) * curvatureNorm > 0.2) {
    return std::nullopt;
  }
  const double halfPathSquared = 0.5 * s * s;
  double ds = 0.;
  double displacementScale = halfPathSquared;
  if constexpr (Quadratic) {
    // Retain the s*ds*k and ds^2*k/2 terms in the local parabola:
    //   a*ds^2 + b*ds + c = 0,
    //   a = n.k/2, b = n.t + s*n.k, c = (s^2/2)*n.k.
    // The rationalized root is continuous with ds = 0 as k -> 0.
    // A reversed slope or negative discriminant invalidates this branch.
    // This solves the seed tangent plane, not the full cylinder or helix.
    const double a = 0.5 * normalCurvature;
    const double b = incidence + s * normalCurvature;
    const double c = halfPathSquared * normalCurvature;
    const double discriminant = b * b - 4. * a * c;
    if (discriminant < 0. || b * incidence <= 0.) {
      return std::nullopt;
    }
    ds = -2. * c / (b + std::copysign(std::sqrt(discriminant), b));
    displacementScale = 0.5 * (s + ds) * (s + ds);
  } else {
    ds = -halfPathSquared * normalCurvature / incidence;
  }
  if (!std::isfinite(ds) || std::abs(ds) > 0.2 * std::abs(s)) {
    return std::nullopt;
  }
  const Vector3 point = seed + displacementScale * curvature + ds * direction;
  if (!point.allFinite()) {
    return std::nullopt;
  }
  return std::pair<Vector3, double>{point, s + ds};
}

}  // namespace Acts::detail
