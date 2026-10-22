// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "OnAxisCylinderIntersection.hpp"

namespace ActsTests {
Acts::MultiIntersection3D intersectOnAxisCylinder(
    const Acts::CylinderSurface& surface, const Acts::GeometryContext& gctx,
    const Acts::Vector3& position, const Acts::Vector3& direction,
    const Acts::BoundaryTolerance& boundaryTolerance, double tolerance) {
  const auto& trf = surface.localToGlobalTransform(gctx);
  const Acts::Vector3 axis = trf.linear().col(2);
  const Acts::Vector3 delta = position - trf.translation();
  const Acts::Vector3 perpendicular = delta.cross(axis);
  const double a = direction.cross(axis).squaredNorm();
  if (!surface.bounds().coversFullAzimuth() ||
      !(perpendicular.array() == 0).all() || a == 0 ||
      (!boundaryTolerance.isNone() && !boundaryTolerance.isInfinite())) {
    return surface.intersect(gctx, position, direction, boundaryTolerance,
                             tolerance);
  }
  const double t =
      surface.bounds().get(Acts::CylinderBounds::eR) / std::sqrt(a);
  // Leave non-finite/underflow extremes to the original solver as well.
  if (!std::isfinite(t) || t == 0) {
    return surface.intersect(gctx, position, direction, boundaryTolerance,
                             tolerance);
  }
  auto intersection = [&](double path) {
    Acts::Vector3 hit = position + path * direction;
    auto status = std::abs(path) < std::abs(tolerance)
                      ? Acts::IntersectionStatus::onSurface
                      : Acts::IntersectionStatus::reachable;
    if (boundaryTolerance.isNone()) {
      const double z = (hit - trf.translation()).dot(axis);
      const double halfZ =
          surface.bounds().get(Acts::CylinderBounds::eHalfLengthZ) + tolerance;
      if (!(std::abs(z) < std::abs(halfZ))) {
        status = Acts::IntersectionStatus::unreachable;
      }
    }
    return Acts::Intersection3D(hit, path, status);
  };
  return Acts::MultiIntersection3D(intersection(-t), intersection(t));
}
}  // namespace ActsTests
