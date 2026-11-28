// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/detail/TwoTrackVertexSeed.hpp"

#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>

namespace Acts::detail {
namespace {
struct Circle {
  Vector3 position;
  Vector2 center;
  double radius;
  double cotTheta;

  // ACTS positive q*B turns clockwise: signed radius = -pT/(q*B).
  std::pair<Vector3, double> at(const Vector2& xy) const {
    const Vector2 start = position.head<2>() - center;
    const Vector2 end = xy - center;
    const double angle =
        std::atan2(start.x() * end.y() - start.y() * end.x(), start.dot(end));
    const double arc = radius * angle;
    return {Vector3(xy.x(), xy.y(), position.z() + arc * cotTheta),
            std::abs(arc)};
  }
};
}  // namespace

Result<TwoTrackVertexSeed> estimateTwoTrackVertexSeed(
    const BoundTrackParameters& first, const BoundTrackParameters& second,
    const GeometryContext& gctx, double fieldZ) {
  const Vector3 p1 = first.position(gctx);
  const Vector3 p2 = second.position(gctx);
  const Vector3 d1 = first.direction();
  const Vector3 d2 = second.direction();
  if (!p1.allFinite() || !p2.allFinite() || !d1.allFinite() ||
      !d2.allFinite() || !std::isfinite(fieldZ)) {
    return VertexingError::InvalidInput;
  }
  if (fieldZ == 0.) {
    const double dot = d1.dot(d2);
    const double determinant = 1. - dot * dot;
    if (determinant < 1e-12) {
      return VertexingError::SeedingError;
    }
    const Vector3 delta = p2 - p1;
    const double s = (delta.dot(d1) - dot * delta.dot(d2)) / determinant;
    const double t = (dot * delta.dot(d1) - delta.dot(d2)) / determinant;
    const Vector3 a = p1 + s * d1;
    const Vector3 b = p2 + t * d2;
    return TwoTrackVertexSeed{(a + b) / 2., (a - b).norm()};
  }

  const double q1 = first.parameters()[eBoundQOverP];
  const double q2 = second.parameters()[eBoundQOverP];
  const double transverse1 = d1.head<2>().norm();
  const double transverse2 = d2.head<2>().norm();
  if (!std::isfinite(q1) || !std::isfinite(q2) || q1 == 0. || q2 == 0. ||
      first.particleHypothesis().absoluteCharge() == 0. ||
      second.particleHypothesis().absoluteCharge() == 0. ||
      transverse1 < 1e-12 || transverse2 < 1e-12) {
    return VertexingError::InvalidInput;
  }
  const double r1 = -transverse1 / (q1 * fieldZ);
  const double r2 = -transverse2 / (q2 * fieldZ);
  const Circle c1{p1,
                  p1.head<2>() + r1 * Vector2(-d1.y(), d1.x()) / transverse1,
                  r1, d1.z() / transverse1};
  const Circle c2{p2,
                  p2.head<2>() + r2 * Vector2(-d2.y(), d2.x()) / transverse2,
                  r2, d2.z() / transverse2};
  if (!c1.center.allFinite() || !c2.center.allFinite()) {
    return VertexingError::InvalidInput;
  }
  const double a = std::abs(r1);
  const double b = std::abs(r2);
  const Vector2 delta = c2.center - c1.center;
  const double distance = delta.norm();
  const double tolerance =
      64. * std::numeric_limits<double>::epsilon() * std::max({a, b, distance});
  if (distance <= tolerance) {
    return VertexingError::SeedingError;
  }
  const Vector2 axis = delta / distance;
  std::array<std::pair<Vector2, Vector2>, 2> candidates;
  std::size_t count = 1;
  if (distance > a + b) {
    candidates[0] = {c1.center + a * axis, c2.center - b * axis};
  } else if (distance < std::abs(a - b)) {
    const double sign = a > b ? 1. : -1.;
    candidates[0] = {c1.center + sign * a * axis, c2.center + sign * b * axis};
  } else {
    const double along =
        ((a - b) * (a + b) + distance * distance) / (2. * distance);
    const double height = std::sqrt(std::max(0., (a - along) * (a + along)));
    const Vector2 base = c1.center + along * axis;
    const Vector2 offset = height * Vector2(-axis.y(), axis.x());
    candidates[0] = {base + offset, base + offset};
    candidates[1] = {base - offset, base - offset};
    count = 2;
  }
  TwoTrackVertexSeed best;
  best.separation = std::numeric_limits<double>::infinity();
  double bestArc = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < count; ++i) {
    const auto [point1, arc1] = c1.at(candidates[i].first);
    const auto [point2, arc2] = c2.at(candidates[i].second);
    const double separation = (point1 - point2).norm();
    const double arc = arc1 + arc2;
    if (separation < best.separation - tolerance ||
        (std::abs(separation - best.separation) <= tolerance &&
         arc < bestArc)) {
      best = {(point1 + point2) / 2., separation};
      bestArc = arc;
    }
  }
  if (!best.position.allFinite() || !std::isfinite(best.separation)) {
    return VertexingError::NumericFailure;
  }
  return best;
}
}  // namespace Acts::detail
