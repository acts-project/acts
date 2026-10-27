// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cmath>

namespace Acts::detail {

/// First-order local-curvature correction of straight surface intersections.
///
/// For each unbounded ray crossing at s, displace by delta = s*s*curvature/2
/// and adjust the path by ds = -normal.dot(delta)/normal.dot(direction).
/// Bounds are checked at the corrected position, so a bounded ray miss can
/// become a candidate. Cylinder crossings are projected back onto the cylinder
/// after the tangent-plane correction. Intersection indices are preserved.
///
/// This is a local approximation, not a helix intersection: it requires a ray
/// crossing seed. Unsupported surface types, grazing incidence, or large
/// corrections fall back to the bounded straight intersection. The limits
/// below bound the perturbation, not its error. Neighbor coverage and policy
/// validity remain the caller's responsibility.
///
/// @param surface Surface to intersect (correction supports plane, disc, cylinder)
/// @param gctx Geometry context
/// @param position Starting global position
/// @param direction Unit propagation direction
/// @param curvature Global geometric curvature, unchanged under path reversal
/// @param boundaryTolerance Bounds tolerance at the corrected crossing
/// @param tolerance On-surface tolerance
/// @return Corrected intersections, with original ray branch indices
inline MultiIntersection3D intersectWithCurvature(
    const Surface& surface, const GeometryContext& gctx,
    const Vector3& position, const Vector3& direction, const Vector3& curvature,
    const BoundaryTolerance& boundaryTolerance,
    double tolerance = s_onSurfaceTolerance) {
  const auto straight = [&] {
    return surface.intersect(gctx, position, direction, boundaryTolerance,
                             tolerance);
  };
  const bool supported = surface.type() == Surface::Plane ||
                         surface.type() == Surface::Disc ||
                         surface.type() == Surface::Cylinder;
  if (curvature.isZero(0.) || !supported || !curvature.allFinite()) {
    return straight();
  }
  const auto seeds = surface.intersect(
      gctx, position, direction, BoundaryTolerance::Infinite(), tolerance);
  const auto correct = [&](IntersectionIndex index) -> Intersection3D {
    const auto& seed = seeds.at(index);
    if (!seed.isValid()) {
      return seed;
    }
    // Arrival must still be decided by the actual position and physical bounds.
    // Do not turn a predicted crossing into an observed surface hit.
    if (seed.status() == IntersectionStatus::onSurface) {
      return straight().at(index);
    }
    const double s = seed.pathLength();
    const Vector3 normal = surface.normal(gctx, seed.position(), direction);
    const double incidence = normal.dot(direction);
    if (std::abs(incidence) < 1e-6 || std::abs(s) * curvature.norm() > 0.2) {
      return straight().at(index);
    }
    const Vector3 delta = (0.5 * s * s) * curvature;
    const double ds = -normal.dot(delta) / incidence;
    if (!std::isfinite(ds) || std::abs(ds) > 0.2 * std::abs(s)) {
      return straight().at(index);
    }
    Vector3 crossing = seed.position() + delta + ds * direction;
    if (surface.type() == Surface::Cylinder) {
      const auto& transform = surface.localToGlobalTransform(gctx);
      Vector3 local = transform.inverse() * crossing;
      const double radius = static_cast<const CylinderBounds&>(surface.bounds())
                                .get(CylinderBounds::eR);
      const double r = local.head<2>().norm();
      if (!(r > 0.) || !std::isfinite(r)) {
        return straight().at(index);
      }
      local.head<2>() *= radius / r;
      crossing = transform * local;
    }
    if (!crossing.allFinite()) {
      return straight().at(index);
    }
    const bool inside = surface.isOnSurface(gctx, crossing, direction,
                                            boundaryTolerance, tolerance);
    return {crossing, s + ds,
            inside ? IntersectionStatus::reachable
                   : IntersectionStatus::unreachable};
  };
  if (seeds.size() == 1) {
    return MultiIntersection3D(correct(0));
  }
  return MultiIntersection3D(correct(0), correct(1));
}

}  // namespace Acts::detail
