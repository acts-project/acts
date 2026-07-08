// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>

#include <Eigen/Eigenvalues>

namespace Acts::detail {

/// Result of a circle fit in a plane.
struct CircleFit {
  /// Circle center.
  Vector2 center = Vector2::Zero();
  /// Circle radius.
  double radius = 0.;
  /// Whether a finite circle could be fitted. `false` signals the straight-line
  /// (infinite-radius) limit, i.e. (near-)collinear input.
  bool valid = false;
};

/// Algebraic circle fit using the Taubin method over the transverse `(x, y)`
/// projection of the given points.
///
/// The Taubin method fits the general conic `A(x^2+y^2) + B*x + C*y + D = 0`
/// with `A` free to reach zero, so it degrades gracefully towards the
/// straight-line (infinite-radius) limit instead of becoming ill-conditioned.
/// The minimizer of the Taubin cost is the eigenvector of the smallest
/// generalized eigenvalue of a small 3x3 problem, assembled from the
/// centroid-relative moments and solved with Eigen. Straightness is detected
/// from the fitted curvature (radius large compared to the point spread), not
/// from any field threshold, so both weakly bent and genuinely straight tracks
/// fall through to `valid = false`.
///
/// Optional per-point weights turn the fit into a weighted algebraic fit: the
/// centroid and all moments become weighted means. Weights act as relative
/// (e.g. inverse-variance) factors on the transverse residuals; an empty span
/// means uniform weights.
///
/// @param points the points whose transverse projection is fitted
/// @param weights optional per-point weights (empty span = uniform)
/// @return the fitted circle (with `valid = false` if effectively straight)
inline CircleFit fitCircleTaubin(std::span<const Vector3> points,
                                 std::span<const double> weights = {}) {
  CircleFit result;
  const std::size_t n = points.size();
  if (n < 3) {
    return result;
  }
  assert((weights.empty() || weights.size() == n) &&
         "weights must be empty or match the number of points");

  // Per-point weight: an empty span means uniform weights of one.
  const auto weightAt = [&](std::size_t i) {
    return weights.empty() ? 1. : weights[i];
  };

  // Weighted centroid of the transverse projection.
  double sumW = 0.;
  double meanX = 0.;
  double meanY = 0.;
  for (std::size_t i = 0; i < n; ++i) {
    const double w = weightAt(i);
    sumW += w;
    meanX += w * points[i].x();
    meanY += w * points[i].y();
  }
  if (!(sumW > 0.)) {
    // Non-positive total weight: no fit can be formed.
    return result;
  }
  const double invW = 1. / sumW;
  meanX *= invW;
  meanY *= invW;

  // Centroid-relative weighted moments (weighted means), with z = x^2 + y^2.
  double mxx = 0., myy = 0., mxy = 0., mxz = 0., myz = 0., mzz = 0., mz = 0.;
  for (std::size_t i = 0; i < n; ++i) {
    const double w = weightAt(i);
    const double x = points[i].x() - meanX;
    const double y = points[i].y() - meanY;
    const double z = x * x + y * y;
    mxx += w * x * x;
    myy += w * y * y;
    mxy += w * x * y;
    mxz += w * x * z;
    myz += w * y * z;
    mzz += w * z * z;
    mz += w * z;
  }
  mxx *= invW;
  myy *= invW;
  mxy *= invW;
  mxz *= invW;
  myz *= invW;
  mzz *= invW;
  mz *= invW;
  if (!(mz > 0.)) {
    // All points coincide.
    return result;
  }

  // Reduced Taubin generalized eigenproblem `M a = lambda N a` for the reduced
  // parameter vector a = (A, B, C) (the offset D is eliminated as D = -A*mz).
  SquareMatrix3 m;
  m << mzz - mz * mz, mxz, myz,  //
      mxz, mxx, mxy,             //
      myz, mxy, myy;
  SquareMatrix3 nMat = SquareMatrix3::Zero();
  nMat(0, 0) = 4. * mz;
  nMat(1, 1) = 1.;
  nMat(2, 2) = 1.;

  Eigen::GeneralizedSelfAdjointEigenSolver<SquareMatrix3> solver(m, nMat);
  if (solver.info() != Eigen::Success) {
    return result;
  }
  // Eigenvalues are sorted ascending, so the first eigenvector minimizes the
  // Taubin cost.
  const Vector3 a = solver.eigenvectors().col(0);
  const double coeffA = a(0);
  const double coeffB = a(1);
  const double coeffC = a(2);
  if (coeffA == 0.) {
    // Perfectly straight: no finite circle.
    return result;
  }

  const double coeffD = -coeffA * mz;
  const double r2 = (coeffB * coeffB + coeffC * coeffC - 4. * coeffA * coeffD) /
                    (4. * coeffA * coeffA);
  if (!std::isfinite(r2) || r2 <= 0.) {
    return result;
  }
  const double radius = std::sqrt(r2);

  // Straightness by curvature: a radius much larger than the point spread
  // (sqrt(mz) is the RMS distance from the centroid) is treated as a line.
  constexpr double maxRadiusToSpread = 1e6;
  if (radius > maxRadiusToSpread * std::sqrt(mz)) {
    return result;
  }

  const Vector2 centerRel(-coeffB / (2. * coeffA), -coeffC / (2. * coeffA));
  result.center = centerRel + Vector2(meanX, meanY);
  result.radius = radius;
  result.valid = std::isfinite(result.radius) && result.radius > 0.;
  return result;
}

/// Refine a circle fit geometrically by minimizing the true radial residuals
/// `sum_i (|p_i - center| - R)^2` with a few Gauss-Newton steps, seeded by an
/// algebraic (e.g. Taubin) fit. Operates on the transverse `(x, y)` projection
/// of the given points and leaves an invalid seed untouched.
///
/// Optional per-point weights scale each point's contribution to the normal
/// equations, matching the weighting of @ref fitCircleTaubin; an empty span
/// means uniform weights.
///
/// @param fit the circle fit to refine in place
/// @param points the points whose transverse projection is fitted
/// @param iterations the maximum number of Gauss-Newton iterations
/// @param weights optional per-point weights (empty span = uniform)
inline void refineCircleGeometric(CircleFit& fit,
                                  std::span<const Vector3> points,
                                  const std::size_t iterations,
                                  std::span<const double> weights = {}) {
  if (!fit.valid) {
    return;
  }
  assert((weights.empty() || weights.size() == points.size()) &&
         "weights must be empty or match the number of points");

  // Per-point weight: an empty span means uniform weights of one.
  const auto weightAt = [&](std::size_t i) {
    return weights.empty() ? 1. : weights[i];
  };

  for (std::size_t it = 0; it < iterations; ++it) {
    SquareMatrix3 jtj = SquareMatrix3::Zero();
    Vector3 jtr = Vector3::Zero();
    for (std::size_t i = 0; i < points.size(); ++i) {
      const Vector2 d = points[i].head<2>() - fit.center;
      const double dist = d.norm();
      if (dist < std::numeric_limits<double>::epsilon()) {
        continue;
      }
      const double w = weightAt(i);
      const double residual = dist - fit.radius;
      // Jacobian of the residual w.r.t. (cx, cy, R).
      const Vector3 j(-d.x() / dist, -d.y() / dist, -1.);
      jtj += w * j * j.transpose();
      jtr += w * j * residual;
    }
    const Vector3 delta = jtj.ldlt().solve(-jtr);
    if (!delta.allFinite()) {
      break;
    }
    fit.center.x() += delta.x();
    fit.center.y() += delta.y();
    fit.radius += delta.z();
    if (!(fit.radius > 0.)) {
      fit.valid = false;
      return;
    }
    if (delta.norm() < 1e-9) {
      break;
    }
  }
}

}  // namespace Acts::detail
