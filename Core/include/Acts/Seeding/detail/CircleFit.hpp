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
#include <optional>
#include <span>

#include <Eigen/Eigenvalues>

namespace Acts::detail {

/// A circle in a plane, as fitted to a set of points.
struct CircleFit {
  /// Circle center.
  Vector2 center = Vector2::Zero();
  /// Circle radius.
  double radius = 0;
};

/// Algebraic circle fit (Taubin) over the transverse `(x, y)` projection.
///
/// Fits `A(x^2+y^2) + B*x + C*y + D = 0` with `A` free to reach zero, so
/// near-collinear input degrades to a line, reported as `std::nullopt`,
/// instead of becoming ill-conditioned. Weights are relative factors on the
/// residuals; an empty span means uniform.
///
/// @param points the points whose transverse projection is fitted
/// @param weights optional per-point weights (empty span = uniform)
/// @return the fitted circle, or `std::nullopt` if no finite circle is defined
inline std::optional<CircleFit> fitCircleTaubin(
    std::span<const Vector3> points, std::span<const double> weights = {}) {
  const std::size_t n = points.size();
  if (n < 3) {
    return std::nullopt;
  }
  assert((weights.empty() || weights.size() == n) &&
         "weights must be empty or match the number of points");

  const auto weightAt = [&](std::size_t i) {
    return weights.empty() ? 1 : weights[i];
  };

  // Weighted centroid of the transverse projection.
  double sumW = 0;
  double meanX = 0;
  double meanY = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double w = weightAt(i);
    sumW += w;
    meanX += w * points[i].x();
    meanY += w * points[i].y();
  }
  if (sumW <= 0) {
    return std::nullopt;
  }
  const double invW = 1 / sumW;
  meanX *= invW;
  meanY *= invW;

  // Centroid-relative weighted moments, z = x^2 + y^2.
  double mxx = 0;
  double myy = 0;
  double mxy = 0;
  double mxz = 0;
  double myz = 0;
  double mzz = 0;
  double mz = 0;
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
  if (mz <= 0) {
    // All points coincide.
    return std::nullopt;
  }

  // Taubin eigenproblem `M a = lambda N a` for a = (A, B, C), with D = -A*mz.
  SquareMatrix3 m;
  m << mzz - mz * mz, mxz, myz,  //
      mxz, mxx, mxy,             //
      myz, mxy, myy;
  SquareMatrix3 nMat = SquareMatrix3::Zero();
  nMat(0, 0) = 4 * mz;
  nMat(1, 1) = 1;
  nMat(2, 2) = 1;

  Eigen::GeneralizedSelfAdjointEigenSolver<SquareMatrix3> solver(m, nMat);
  if (solver.info() != Eigen::Success) {
    return std::nullopt;
  }
  // Ascending eigenvalues, so the first eigenvector minimizes the Taubin cost.
  const Vector3 a = solver.eigenvectors().col(0);
  const double coeffA = a(0);
  const double coeffB = a(1);
  const double coeffC = a(2);
  if (coeffA == 0) {
    // Perfectly straight: no finite circle.
    return std::nullopt;
  }

  const double coeffD = -coeffA * mz;
  const double r2 = (coeffB * coeffB + coeffC * coeffC - 4 * coeffA * coeffD) /
                    (4 * coeffA * coeffA);
  if (!std::isfinite(r2) || r2 <= 0) {
    return std::nullopt;
  }
  const double radius = std::sqrt(r2);

  // A radius far beyond the point spread sqrt(mz) is a line.
  if (constexpr double maxRadiusToSpread = 1e6;
      radius > maxRadiusToSpread * std::sqrt(mz)) {
    return std::nullopt;
  }

  const Vector2 centerRel(-coeffB / (2 * coeffA), -coeffC / (2 * coeffA));
  return CircleFit{centerRel + Vector2(meanX, meanY), radius};
}

/// Refine a circle fit by minimizing the radial residuals
/// `sum_i (|p_i - center| - R)^2` with Gauss-Newton steps, on the transverse
/// `(x, y)` projection. Weights as in @ref fitCircleTaubin.
///
/// @param fit the circle fit to refine, e.g. from @ref fitCircleTaubin
/// @param points the points whose transverse projection is fitted
/// @param iterations the maximum number of Gauss-Newton iterations
/// @param weights optional per-point weights (empty span = uniform)
/// @return the refined circle, or `std::nullopt` if a step drives the radius
///         non-positive, i.e. the refinement collapses the circle
inline std::optional<CircleFit> refineCircleGeometric(
    CircleFit fit, std::span<const Vector3> points,
    const std::size_t iterations, std::span<const double> weights = {}) {
  assert((weights.empty() || weights.size() == points.size()) &&
         "weights must be empty or match the number of points");

  const auto weightAt = [&](std::size_t i) {
    return weights.empty() ? 1 : weights[i];
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
      const Vector3 j(-d.x() / dist, -d.y() / dist, -1);
      jtj += w * j * j.transpose();
      jtr += w * j * residual;
    }
    const Vector3 delta = jtj.ldlt().solve(-jtr);
    if (!delta.allFinite()) {
      // Singular normal equations: keep the fit as it is.
      break;
    }
    fit.center.x() += delta.x();
    fit.center.y() += delta.y();
    fit.radius += delta.z();
    if (fit.radius <= 0) {
      return std::nullopt;
    }
    if (delta.norm() < 1e-9) {
      break;
    }
  }
  return fit;
}

}  // namespace Acts::detail
