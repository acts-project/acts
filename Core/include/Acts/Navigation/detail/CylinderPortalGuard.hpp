// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/detail/CylinderPortalBounds.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Acts::detail {
// Experimental: homogeneous axial field, fixed kinematics and immutable policy
// geometry. A local curvature sample is not a bound on a general future field.
// Supporting-plane waypoints become Portal hits only after finite owner-face
// checks. Sensor selection and native ray intersection indices are preserved.
struct CylinderPortalGuardBase : CylinderPortalBounds {
  bool radialOverride = false;
  bool positiveStep = false;
  const Portal* pendingPortal = nullptr;
  double pendingFlight = std::numeric_limits<double>::infinity();

  template <typename NavigationState>
  void beforeStep(const GeometryContext& /*gctx*/, const NavigationState& nav,
                  const Vector3& position, const Vector3& direction,
                  double proposedStep) {
    pendingPortal = nullptr;
    pendingFlight = std::numeric_limits<double>::infinity();
    // Negative steps converge to an already selected native intersection root.
    if (!(proposedStep > 0.) || cachedVolume != nav.currentVolume ||
        !nav.stream.isValid()) {
      return;
    }
    if (nestedContainer) {
      // Each child endcap is a real interface too. Its exact axial flight also
      // expires any radial guard that was inapplicable outside the child's z
      // interval. Finite Portal bounds still decide arrival below.
      if (direction.z() == 0.) {
        return;
      }
      for (std::size_t i = 0; i < containerDiscCount; ++i) {
        const auto [portal, z] = containerDiscs[i];
        const double flight = (z - position.z()) / direction.z();
        if (flight > nav.options.surfaceTolerance && flight <= proposedStep &&
            flight < pendingFlight) {
          pendingPortal = portal;
          pendingFlight = flight;
        }
      }
      return;
    }
    const bool forward = direction.z() > 0.;
    const Portal* face = forward ? upperDisc : lowerDisc;
    if (face == nullptr) {
      return;
    }
    const double speed = std::abs(direction.z());
    const double distance =
        forward ? upperZ - position.z() : position.z() - lowerZ;
    // The common non-crossing path avoids both division and surface calls.
    if (distance > nav.options.surfaceTolerance * speed &&
        distance <= proposedStep * speed) {
      pendingPortal = face;
      pendingFlight = distance / speed;
    }
  }

  template <typename NavigationState>
  NavigationTarget radialAlternative(const GeometryContext& gctx,
                                     const NavigationState& nav,
                                     const Vector3& position,
                                     const Vector3& direction,
                                     const NavigationTarget& active,
                                     double proposedStep,
                                     double activeFlight) const {
    // Limit this prototype to portal competition. Sensor lookup, selection,
    // bounds and arrival remain exactly the original implementation.
    if (!(proposedStep > 0.) ||
        !(activeFlight > nav.options.surfaceTolerance) ||
        !active.isPortalTarget() || cachedVolume != nav.currentVolume ||
        innerCylinder == nullptr || outerCylinder == nullptr) {
      return NavigationTarget::None();
    }
    const Vector2 p = position.head<2>() - center;
    const Vector2 d = direction.head<2>();
    const double radialDot = p.dot(d);
    const Portal* alternative = radialDot < 0. ? innerCylinder : outerCylinder;
    if (&active.portal() == alternative) {
      return NavigationTarget::None();
    }
    const double r2 = p.squaredNorm();
    if (radialDot < 0.) {
      // No root/sqrt needed until an inward tangent can actually reach the
      // inner cylinder. This detects a root that was absent at volume entry.
      const double cross = p.x() * d.y() - p.y() * d.x();
      if (!(r2 > innerRadius2) ||
          cross * cross >= innerRadius2 * d.squaredNorm()) {
        return NavigationTarget::None();
      }
    } else if (!(r2 < outerRadius2)) {
      return NavigationTarget::None();
    }
    const auto intersections = alternative->surface().intersect(
        gctx, position, direction, BoundaryTolerance::None(),
        nav.options.surfaceTolerance);
    const auto [intersection, index] = intersections.closestForwardWithIndex();
    if (!intersection.isValid() ||
        intersection.pathLength() <= nav.options.surfaceTolerance ||
        intersection.pathLength() >= activeFlight) {
      return NavigationTarget::None();
    }
    return NavigationTarget(intersection, index, *alternative,
                            BoundaryTolerance::None());
  }

  template <typename NavigationState>
  double radialStepLimit(const NavigationState& nav, const Vector3& position,
                         const Vector3& direction, const Vector3& curvature,
                         const NavigationTarget& active,
                         double proposedStep) const {
    if (!(proposedStep > 0.) || cachedVolume != nav.currentVolume ||
        innerCylinder == nullptr || outerCylinder == nullptr) {
      return proposedStep;
    }
    const Vector2 p = position.head<2>() - center;
    const Vector2 d = direction.head<2>();
    const double c = p.squaredNorm();
    const double b = p.dot(d);
    const double a = d.squaredNorm();
    const double radius = std::sqrt(outerRadius2);
    const double tolerance = 2. * radius * nav.options.surfaceTolerance;
    if (c < innerRadius2 - tolerance || c > outerRadius2 + tolerance) {
      return proposedStep;  // Already outside owner: not a containment
                            // certificate.
    }
    const double k2 = curvature.squaredNorm();
    const double l2 = proposedStep * proposedStep;
    const double bend2 = outerRadius2 * k2 * l2 * l2;
    double cap = proposedStep;
    const bool targetInner =
        active.isPortalTarget() && &active.portal() == innerCylinder;
    const bool targetOuter =
        active.isPortalTarget() && &active.portal() == outerCylinder;
    const double innerClearance = c + 2. * b * proposedStep - innerRadius2;
    if (!targetInner &&
        !(innerClearance > 0. && innerClearance * innerClearance > bend2)) {
      const double C = std::max(0., c - innerRadius2);
      const double alpha = radius * std::sqrt(k2);
      double h = std::numeric_limits<double>::infinity();
      if (alpha == 0.) {
        if (b < 0.) {
          h = C / (-2. * b);
        }
      } else {
        const double root = std::sqrt(b * b + alpha * C);
        // At an entering inner face C=0,b>0, keep the nonzero return root.
        h = b >= 0. ? (b + root) / alpha : C / (root - b);
      }
      if (nestedContainer) {
        // The supporting cylinder does not exist beyond the child's finite
        // z extent. A real radial crossing cannot occur before BOTH the radial
        // certificate and entry into that axial interval. This avoids an
        // asymptotic approach to a nonexistent supporting-cylinder crossing.
        double axialEntry = 0.;
        if (position.z() < innerMinZ) {
          axialEntry = direction.z() > 0.
                           ? (innerMinZ - position.z()) / direction.z()
                           : std::numeric_limits<double>::infinity();
        } else if (position.z() > innerMaxZ) {
          axialEntry = direction.z() < 0.
                           ? (innerMaxZ - position.z()) / direction.z()
                           : std::numeric_limits<double>::infinity();
        }
        h = std::max(h, axialEntry);
      }
      if (h > 0.) {
        cap = std::min(cap, h);
      }
    }
    // In the homogeneous axial field a=|direction_xy|^2 is constant.
    const double outerClearance =
        outerRadius2 - c - 2. * b * proposedStep - a * l2;
    if (!targetOuter &&
        !(outerClearance > 0. && outerClearance * outerClearance > bend2)) {
      const double C = std::max(0., outerRadius2 - c);
      const double A = a + radius * std::sqrt(k2);
      if (A > 0.) {
        const double root = std::sqrt(b * b + A * C);
        const double h = b >= 0. ? C / (root + b) : (root - b) / A;
        if (h > 0.) {
          cap = std::min(cap, h);
        }
      }
    }
    // A zero bound at exact tangency is unresolved by these inequalities;
    // no artificial minimum step is advertised as conservative. The physical
    // face-arrival hook handles outgoing on-face states after positive steps.
    return cap;
  }

  template <typename NavigationState>
  const Portal* physicalRadialExit(const GeometryContext& gctx,
                                   const NavigationState& nav,
                                   const Vector3& position,
                                   const Vector3& direction,
                                   bool steppedForward) const {
    if (!steppedForward || cachedVolume != nav.currentVolume ||
        innerCylinder == nullptr || outerCylinder == nullptr) {
      return nullptr;
    }
    const Vector2 p = position.head<2>() - center;
    const double c = p.squaredNorm();
    const double b = p.dot(direction.head<2>());
    if (b == 0.) {
      return nullptr;
    }
    const Portal* portal = b < 0. ? innerCylinder : outerCylinder;
    const double r2 = b < 0. ? innerRadius2 : outerRadius2;
    const double radius = std::sqrt(r2);
    const double epsilon = nav.options.surfaceTolerance;
    if (std::abs(c - r2) > 2. * radius * epsilon + epsilon * epsilon ||
        !cachedVolume->inside(gctx, position, epsilon) ||
        !portal->surface().isOnSurface(gctx, position, direction,
                                       BoundaryTolerance::None(), epsilon)) {
      return nullptr;
    }
    return portal;
  }

  template <typename NavigationState>
  bool isPhysicalExit(const GeometryContext& gctx, const NavigationState& nav,
                      const Vector3& position, const Vector3& direction,
                      double performedStep) const {
    if (pendingPortal == nullptr ||
        std::abs(performedStep) + nav.options.surfaceTolerance <
            pendingFlight) {
      return false;
    }
    // Shared portal bounds can exceed the actual owner's radial face extent.
    // Do not resolve a support waypoint that is outside the owner's real face.
    return nav.currentVolume == cachedVolume &&
           cachedVolume->inside(gctx, position, nav.options.surfaceTolerance) &&
           pendingPortal->surface().isOnSurface(gctx, position, direction,
                                                BoundaryTolerance::None(),
                                                nav.options.surfaceTolerance);
  }
};

struct CylinderPortalGuard : CylinderPortalGuardBase {
  const Portal* certificateTarget = nullptr;
  double certificateStart = -std::numeric_limits<double>::infinity();
  double certificateEnd = -std::numeric_limits<double>::infinity();
  bool certificateValid = false;
  double innerEnd = 0., outerEnd = 0., discEnd = 0.;
  const Portal* forwardDisc = nullptr;
  std::array<double, 4> containerDiscEnds{};
  double surfaceTolerance = 0.;
  double outerRadius = 0.;

  template <typename NavigationState>
  bool needsCheck(const GeometryContext& /*gctx*/, const NavigationState& nav,
                  const Vector3& position, const Vector3& direction,
                  const Vector3& curvature, const NavigationTarget& active,
                  const ConstrainedStep& stepSize, double path) {
    pendingPortal = nullptr;
    pendingFlight = std::numeric_limits<double>::infinity();
    const auto* descriptor = CylinderPortalBounds::get(nav.stream);
    if (nav.currentVolume == nullptr || !nav.stream.isValid() ||
        descriptor == nullptr ||
        descriptor->cachedVolume != nav.currentVolume) {
      certificateValid = false;
      return false;
    }
    const Portal* target = active.isPortalTarget() ? &active.portal() : nullptr;
    if (cachedVolume != nav.currentVolume) {
      load(*descriptor, nav.options.surfaceTolerance);
    }
    if (!certificateValid || path < certificateStart) {
      refresh(position, direction, curvature, path);
      certificateTarget = target;
      select(target);
    } else if (target != certificateTarget) {
      certificateTarget = target;
      select(target);
    }
    // The user constraint is a cheap upper bound on any actual positive step.
    // Avoid the full min/sign calculation on the common, safely interior path.
    if (path + std::abs(stepSize.value(ConstrainedStep::Type::User)) <
        certificateEnd) {
      return false;
    }
    const double proposed = stepSize.value();
    if (!(proposed > 0.) || path + proposed < certificateEnd) {
      return false;
    }
    // Recompute from the actual state when the old conservative bound expires.
    // This is often still safe: the bound is a lower bound, not a crossing.
    refresh(position, direction, curvature, path);
    select(target);
    if (path + proposed < certificateEnd) {
      return false;
    }
    // The original guard will cap/retarget. Refresh after that changed flight.
    certificateValid = false;
    return true;
  }

  void select(const Portal* target) {
    constexpr double infinity = std::numeric_limits<double>::infinity();
    certificateEnd = std::min(target == innerCylinder ? infinity : innerEnd,
                              target == outerCylinder ? infinity : outerEnd);
    if (nestedContainer) {
      for (std::size_t i = 0; i < containerDiscCount; ++i) {
        if (target != containerDiscs[i].first) {
          certificateEnd = std::min(certificateEnd, containerDiscEnds[i]);
        }
      }
    } else if (target != forwardDisc) {
      certificateEnd = std::min(certificateEnd, discEnd);
    }
  }

  // Keep rare geometry changes and deadline algebra out of the per-step hot
  // path while retaining their definitions in the propagation translation unit.
  [[gnu::noinline]] void load(const CylinderPortalBounds& bounds,
                              double tolerance) {
    static_cast<CylinderPortalBounds&>(*this) = bounds;
    surfaceTolerance = tolerance;
    outerRadius = std::sqrt(outerRadius2);
    radialOverride = false;
    certificateValid = false;
  }

  [[gnu::noinline]] void refresh(const Vector3& position,
                                 const Vector3& direction,
                                 const Vector3& curvature, double path) {
    certificateStart = path;
    certificateValid = true;
    innerEnd = outerEnd = discEnd = std::numeric_limits<double>::infinity();
    containerDiscEnds.fill(std::numeric_limits<double>::infinity());
    const double k2 = curvature.squaredNorm();
    if (k2 == 0.) {
      return;
    }
    const bool forward = direction.z() > 0.;
    forwardDisc = forward ? upperDisc : lowerDisc;
    if (nestedContainer && direction.z() != 0.) {
      for (std::size_t i = 0; i < containerDiscCount; ++i) {
        const double flight =
            (containerDiscs[i].second - position.z()) / direction.z();
        if (flight > surfaceTolerance) {
          containerDiscEnds[i] = path + flight;
        }
      }
    } else if (forwardDisc != nullptr && direction.z() != 0.) {
      const double distance =
          forward ? upperZ - position.z() : position.z() - lowerZ;
      discEnd = path + std::max(0., distance / std::abs(direction.z()));
    }
    if (innerCylinder == nullptr || outerCylinder == nullptr) {
      return;
    }
    const Vector2 p = position.head<2>() - center;
    const Vector2 d = direction.head<2>();
    const double c = p.squaredNorm();
    const double b = p.dot(d);
    const double a = d.squaredNorm();
    // Write g(s)=|p_xy(s)|^2, c=g(0), b=p_xy.d_xy, a=|d_xy|^2.
    // In a homogeneous axial field, a and K=|curvature| are constant.
    // Before the first owner exit, |p_xy| <= R, hence
    //   -2 R K <= g''(s) = 2a + 2 p_xy.curvature_xy <= 2(a + R K).
    // Integrating twice gives c+2bs-RKs^2 <= g(s) <= c+2bs+(a+RK)s^2.
    // The first roots at the omitted inner/outer radius bound how far we can
    // travel without reconsidering that face. They are refresh deadlines, not
    // predicted intersections. Keep the nonzero inward-return root on entry.
    const double alpha = outerRadius * std::sqrt(k2);
    {
      const double C = std::max(0., c - innerRadius2);
      const double root = std::sqrt(b * b + alpha * C);
      innerEnd = path + (b >= 0. ? (b + root) / alpha : C / (root - b));
      if (nestedContainer) {
        // A real inward child-cylinder crossing requires both radial arrival
        // and entry into the finite axial extent. Child endcaps independently
        // expire the cached interval, including positions inside its radius.
        double axialEntry = 0.;
        if (position.z() < innerMinZ) {
          axialEntry = direction.z() > 0.
                           ? (innerMinZ - position.z()) / direction.z()
                           : std::numeric_limits<double>::infinity();
        } else if (position.z() > innerMaxZ) {
          axialEntry = direction.z() < 0.
                           ? (innerMaxZ - position.z()) / direction.z()
                           : std::numeric_limits<double>::infinity();
        }
        // Native arrival may be just inside the radius within tolerance.
        // Retain the nonzero return bound after an outward entry there.
        innerEnd = c < innerRadius2 - 2. * outerRadius * surfaceTolerance
                       ? std::numeric_limits<double>::infinity()
                       : std::max(innerEnd, path + axialEntry);
      }
    }
    {
      const double C = std::max(0., outerRadius2 - c);
      const double A = a + alpha;
      const double root = std::sqrt(b * b + A * C);
      outerEnd = path + (b >= 0. ? C / (root + b) : (root - b) / A);
    }
  }
};

}  // namespace Acts::detail
