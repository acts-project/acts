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
  // Use unbounded ray seeds: a ray outside the sensor bounds can bend inside.
  // Bounds are applied only after correcting the crossing below.
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
    // Near tangency the root shift is singular, so a perturbation about the
    // ray root is unreliable. |s|*|curvature| estimates the local tangent
    // change; its limit and the relative path-shift limit below are heuristic
    // guards, not bounds on the intersection error.
    if (std::abs(incidence) < 1e-6 || std::abs(s) * curvature.norm() > 0.2) {
      return straight().at(index);
    }
    // Local trajectory model, with unit tangent t and k = d^2 x / ds^2:
    //
    //   x(u) ~= x0 + u*t + (u^2/2)*k.
    //
    // k has units of inverse length. We freeze it at the starting point; this
    // omits its evolution along the trajectory and is not an exact helix model.
    // Let y = x0 + s*t be the ray seed and n its surface normal. Write the
    // corrected path as s + ds. Substitution into the local model gives
    //
    //   x(s + ds) = y + (s^2/2)*k + ds*t + s*ds*k + (ds^2/2)*k.
    //
    // For a regular (non-grazing) root, ds is first order in k. Thus s*ds*k
    // and (ds^2/2)*k are second and third order in k, respectively, and are
    // dropped. With Delta = (s^2/2)*k, the linearized surface equation is
    //
    //   n . (Delta + ds*t) = 0  =>  ds = -(n . Delta)/(n . t).
    //
    // We retain this first-order path shift: the corrected point is
    // y + Delta + ds*t, and its path estimate is s + ds. Although the local
    // trajectory is quadratic in path length, the intersection correction is
    // only first order in curvature; no quadratic root equation is solved.
    const Vector3 delta = (0.5 * s * s) * curvature;
    const double ds = -normal.dot(delta) / incidence;
    if (!std::isfinite(ds) || std::abs(ds) > 0.2 * std::abs(s)) {
      return straight().at(index);
    }
    Vector3 crossing = seed.position() + delta + ds * direction;
    if (surface.type() == Surface::Cylinder) {
      // The correction lies on the tangent plane, which is the exact support
      // surface for planes/discs but only a local approximation for cylinders.
      // Restore the cylinder radius before checking bounds. Near the seed this
      // projection changes the point at second order in the displacement; it
      // does not make the trajectory or path estimate second-order accurate.
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
