// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/CylinderNavigationPolicy.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <algorithm>
#include <numbers>
#include <sstream>

namespace Acts {

namespace {
std::string_view faceName(CylinderVolumeBounds::Face face) {
  using enum CylinderVolumeBounds::Face;
  switch (face) {
    case PositiveDisc:
      return "PositiveDisc";
    case NegativeDisc:
      return "NegativeDisc";
    case OuterCylinder:
      return "OuterCylinder";
    case InnerCylinder:
      return "InnerCylinder";
    default:
      return "Unknown";
  }
}
}  // namespace

bool CylinderNavigationPolicy::isApplicable(const TrackingVolume& volume,
                                            const Logger& logger) {
  const auto& volumeBounds = volume.volumeBounds();

  if (volumeBounds.type() != VolumeBounds::eCylinder) {
    ACTS_DEBUG("CylinderNavigationPolicy can only be used with "
               << "cylinder volumes, which " << volume.volumeName()
               << " is not");
    return false;
  }

  const auto& bounds =
      dynamic_cast<const CylinderVolumeBounds&>(volume.volumeBounds());

  if (bounds.get(CylinderVolumeBounds::eMinR) == 0) {
    ACTS_DEBUG("CylinderNavigationPolicy can only be used with "
               << "non-zero inner radius, which " << volume.volumeName()
               << " has not");
    return false;
  }

  if (bounds.get(CylinderVolumeBounds::eHalfPhiSector) < std::numbers::pi) {
    ACTS_DEBUG("CylinderNavigationPolicy can only be used with "
               << "full-phi volumes, which " << volume.volumeName()
               << " is not");
    return false;
  }

  // A volume confining sub-volumes also carries their portals, so this
  // implicitly rejects those as well.
  if (volume.portals().size() != 4) {
    ACTS_DEBUG("CylinderNavigationPolicy can only be used with "
               << "volumes with 4 portals, but " << volume.volumeName()
               << " has " << volume.portals().size());
    return false;
  }

  return true;
}

CylinderNavigationPolicy::CylinderNavigationPolicy(const GeometryContext& gctx,
                                                   const TrackingVolume& volume,
                                                   const Logger& logger)
    : m_volume(&volume) {
  ACTS_VERBOSE("CylinderNavigationPolicy constructor for volume "
               << volume.volumeName());
  using enum CylinderVolumeBounds::Face;

  if (!isApplicable(volume, logger)) {
    ACTS_ERROR("CylinderNavigationPolicy is not applicable to volume "
               << volume.volumeName());
    throw std::invalid_argument(
        "CylinderNavigationPolicy is not applicable to volume " +
        volume.volumeName());
  }

  const auto& bounds =
      dynamic_cast<const CylinderVolumeBounds&>(m_volume->volumeBounds());

  double rMin = bounds.get(CylinderVolumeBounds::eMinR);
  m_rMin = rMin;
  m_rMin2 = rMin * rMin;
  double rMax = bounds.get(CylinderVolumeBounds::eMaxR);
  m_rMax = rMax;
  m_rMax2 = rMax * rMax;
  m_alignable = volume.isAlignable();
  m_center = volume.center(gctx);
  m_halfLengthZ = bounds.get(CylinderVolumeBounds::eHalfLengthZ);

  ACTS_VERBOSE("CylinderNavigationPolicy created for volume "
               << volume.volumeName());

  if (!volume.localToGlobalTransform(gctx).linear().isApprox(
          SquareMatrix3::Identity())) {
    m_itransform = volume.globalToLocalTransform(gctx);
  }

  m_shifted = !m_alignable && !m_itransform.has_value();

  // Since the volume does not store the shell assignment, we have to recover
  // this from the raw portals
  for (const auto& portal : m_volume->portals()) {
    if (const auto* cylBounds =
            dynamic_cast<const CylinderBounds*>(&portal.surface().bounds());
        cylBounds != nullptr) {
      if (std::abs(cylBounds->get(CylinderBounds::eR) - rMin) <
          s_onSurfaceTolerance) {
        m_portals.at(toUnderlying(InnerCylinder)) = &portal;
        continue;
      }

      if (std::abs(cylBounds->get(CylinderBounds::eR) - rMax) <
          s_onSurfaceTolerance) {
        m_portals.at(toUnderlying(OuterCylinder)) = &portal;
        continue;
      }
    }

    if (const auto* radBounds =
            dynamic_cast<const RadialBounds*>(&portal.surface().bounds());
        radBounds != nullptr) {
      Transform3 localTransform = m_volume->globalToLocalTransform(gctx) *
                                  portal.surface().localToGlobalTransform(gctx);
      Vector3 localPosition = localTransform.translation();
      double localZ = localPosition.z();
      if (std::abs(localZ - m_halfLengthZ) < s_onSurfaceTolerance) {
        m_portals.at(toUnderlying(PositiveDisc)) = &portal;
        continue;
      }

      if (std::abs(localZ + m_halfLengthZ) < s_onSurfaceTolerance) {
        m_portals.at(toUnderlying(NegativeDisc)) = &portal;
        continue;
      }
    }
  }

  // Validate that the geometric recovery binned a portal onto every face. If
  // not, dump a detailed comparison of every raw portal against the volume's
  // reference geometry so we can see *why* a face was left unassigned (portal
  // dropped because its bounds type / radius / localZ did not match, or two
  // portals collided on the same face).
  const bool complete = std::ranges::all_of(
      m_portals, [](const Portal* portal) { return portal != nullptr; });

  if (!complete) {
    std::ostringstream diag;
    diag << "CylinderNavigationPolicy failed to recover a portal for every "
            "face of volume "
         << volume.volumeName() << "\n";
    diag << "  reference bounds: rMin=" << rMin << " rMax=" << rMax
         << " halfLengthZ=" << m_halfLengthZ
         << " (s_onSurfaceTolerance=" << s_onSurfaceTolerance << ")\n";
    diag << "  volume center (global): " << volume.center(gctx).transpose()
         << "\n";
    diag << "  localToGlobal translation: "
         << volume.localToGlobalTransform(gctx).translation().transpose()
         << "\n";

    diag << "  " << m_volume->portals().size() << " raw portal(s):\n";
    for (const auto& portal : m_volume->portals()) {
      const auto& surface = portal.surface();
      diag << "    - " << surface.toStream(gctx) << "\n";
      if (const auto* cylBounds =
              dynamic_cast<const CylinderBounds*>(&surface.bounds());
          cylBounds != nullptr) {
        const double r = cylBounds->get(CylinderBounds::eR);
        diag << "      CylinderBounds: r=" << r
             << " |r-rMin|=" << std::abs(r - rMin)
             << " |r-rMax|=" << std::abs(r - rMax) << "\n";
      } else if (const auto* radBounds =
                     dynamic_cast<const RadialBounds*>(&surface.bounds());
                 radBounds != nullptr) {
        const Transform3 localTransform =
            m_volume->globalToLocalTransform(gctx) *
            surface.localToGlobalTransform(gctx);
        const double localZ = localTransform.translation().z();
        diag << "      RadialBounds: rMin="
             << radBounds->get(RadialBounds::eMinR)
             << " rMax=" << radBounds->get(RadialBounds::eMaxR)
             << " localZ=" << localZ
             << " |localZ-halfLengthZ|=" << std::abs(localZ - m_halfLengthZ)
             << " |localZ+halfLengthZ|=" << std::abs(localZ + m_halfLengthZ)
             << "\n";
      } else {
        diag << "      (bounds type "
             << static_cast<int>(surface.bounds().type())
             << " is not handled by the recovery logic)\n";
      }
    }

    diag << "  resulting face assignment:\n";
    for (const auto& [i, portal] : enumerate(m_portals)) {
      const auto face = static_cast<CylinderVolumeBounds::Face>(i);
      diag << "    " << faceName(face) << " -> ";
      if (portal == nullptr) {
        diag << "nullptr";
      } else {
        diag << portal->surface().toStream(gctx);
      }
      diag << "\n";
    }

    ACTS_ERROR(diag.str());
    throw std::invalid_argument(diag.str());
  }

  ACTS_VERBOSE("Portal assignment:");
  for (const auto& [i, portal] : enumerate(m_portals)) {
    auto face = static_cast<CylinderVolumeBounds::Face>(i);
    ACTS_VERBOSE("  " << faceName(face) << " -> "
                      << portal->surface().toStream(gctx));
  }
}

void CylinderNavigationPolicy::initializeCandidates(
    [[maybe_unused]] const GeometryContext& gctx,
    const NavigationArguments& args, NavigationPolicyState& /*state*/,
    AppendOnlyNavigationStream& stream, const Logger& logger) const {
  using enum CylinderVolumeBounds::Face;

  ACTS_VERBOSE("CylinderNavigationPolicy::initializeCandidates for volume "
               << m_volume->volumeName()
               << " with gpos: " << args.position.transpose()
               << " gdir: " << args.direction.transpose());

  Vector3 pos;
  Vector3 dir;
  if (m_itransform.has_value()) {
    pos = *m_itransform * args.position;
    dir = m_itransform->linear() * args.direction;
  } else {
    pos = args.position - m_volume->center(gctx);
    dir = args.direction;
  }

  ACTS_VERBOSE("-> lpos: " << pos.transpose() << " ldir: " << dir.transpose());

  auto add = [&](auto face) {
    ACTS_VERBOSE("~~> Adding portal candidate " << faceName(face));
    stream.addPortalCandidate(*m_portals.at(toUnderlying(face)));
  };

  bool hitDisk = false;
  Vector3 diskIntersection;
  double diskDistance3D = std::numeric_limits<double>::max();
  double zDisk = (dir[2] > 0) ? m_halfLengthZ : -m_halfLengthZ;
  if (std::abs(dir[2]) > s_onSurfaceTolerance) {
    ACTS_VERBOSE(
        "-> Not parallel to the disc, see if we're inside the disc radii");

    double t = (zDisk - pos[2]) / dir[2];

    diskIntersection = pos + t * dir;
    double r2 = diskIntersection[0] * diskIntersection[0] +
                diskIntersection[1] * diskIntersection[1];
    if (r2 < m_rMax2 && r2 > m_rMin2 &&
        t > 0) {  // Only consider forward intersections
      hitDisk = true;
      diskDistance3D = t;  // Parameter t is the distance along the ray
    }
  } else {
    ACTS_VERBOSE("-> Parallel to the disc, see if we're inside the disc radii");
  }

  double rpos2 = pos[0] * pos[0] + pos[1] * pos[1];
  ACTS_VERBOSE("-> rpos: " << std::sqrt(rpos2));

  bool hitInner = false;
  if (std::abs(rpos2 - m_rMin2) > s_onSurfaceTolerance) {
    // Find point of closest approach to the origin

    Vector2 dir2 = dir.head<2>().normalized();
    double d = -1 * pos.head<2>().dot(dir2);
    ACTS_VERBOSE("-> d: " << d);

    if (d > 0) {  // Point of closest approach is in the direction of the ray

      double clippedD = d;
      if (hitDisk) {
        // Clip to distance of disk intersection
        double diskIntersectionDistanceXy =
            (diskIntersection.head<2>() - pos.head<2>()).norm();
        clippedD = std::min(d, diskIntersectionDistanceXy);
      }

      Vector2 poc = pos.head<2>() + clippedD * dir2;
      double r2 = poc.dot(poc);
      hitInner = r2 < m_rMin2;
      if (hitInner) {
        add(InnerCylinder);
        // If we hit the inner cylinder before reaching the disk intersection,
        // we can discard the disk as we'll never reach it
        if (hitDisk) {
          // Calculate 3D distance to inner cylinder intersection
          double innerDistance3D =
              d / dir.head<2>().norm();  // Convert 2D distance to 3D parameter
          if (innerDistance3D < diskDistance3D) {
            hitDisk = false;
          }
        }
      }
    }
  }

  // Add disk candidate if we determined it's reachable (not blocked by inner
  // cylinder)
  if (hitDisk) {
    add(dir[2] > 0 ? PositiveDisc : NegativeDisc);
  }

  if (!hitInner && !hitDisk) {
    add(OuterCylinder);
  }
}

void CylinderNavigationPolicy::connect(NavigationDelegate& delegate) const {
  connectDefault<CylinderNavigationPolicy>(delegate);
}

Vector3 CylinderNavigationPolicy::toLocal(const GeometryContext& gctx,
                                          const Vector3& global) const {
  if (m_shifted) {
    return global - m_center;
  }
  if (m_alignable) {
    return m_volume->globalToLocalTransform(gctx) * global;
  }
  return *m_itransform * global;
}

Vector3 CylinderNavigationPolicy::toLocalDirection(
    const GeometryContext& gctx, const Vector3& global) const {
  if (m_shifted) {
    return global;
  }
  if (m_alignable) {
    return m_volume->globalToLocalTransform(gctx).linear() * global;
  }
  return m_itransform->linear() * global;
}

bool CylinderNavigationPolicy::onFace(const Vector3& local,
                                      double tolerance) const {
  const double r2 = local.x() * local.x() + local.y() * local.y();
  const double rMin = std::max(0., m_rMin - tolerance);
  const double rMax = m_rMax + tolerance;
  return r2 >= rMin * rMin && r2 <= rMax * rMax &&
         std::abs(local.z()) <= m_halfLengthZ + tolerance;
}

bool CylinderNavigationPolicy::checkVolumeExit(const GeometryContext& gctx,
                                               const VolumeExitArguments& args,
                                               NavigationTarget& exit,
                                               const Logger& logger) const {
  const Vector3 local = toLocal(gctx, args.position);
  const bool inside = onFace(local, args.tolerance);
  if (inside && (m_rMin == 0. || args.previousPosition == nullptr)) {
    return false;
  }

  // A step that passed the closest approach to the axis can have dipped into
  // the inner cylinder and come back out, also before ending up beyond another
  // face. Along the chord c of the step, that is the case if 0 < p.c < c.c in
  // the transverse plane.
  if (m_rMin > 0. && args.previousPosition != nullptr) {
    const Vector3 chord =
        toLocalDirection(gctx, args.position - *args.previousPosition);
    const double along = local.x() * chord.x() + local.y() * chord.y();
    if (along > 0. && along < chord.x() * chord.x() + chord.y() * chord.y()) {
      if (NavigationTarget crossing = findExit(gctx, args, true);
          !crossing.isNone()) {
        ACTS_VERBOSE("Dipped out of " << m_volume->volumeName());
        exit = crossing;
        return true;
      }
    }
  }
  if (inside) {
    return false;
  }

  NavigationTarget crossing = findExit(gctx, args, false);
  ACTS_VERBOSE("Left " << m_volume->volumeName() << ": "
                       << (crossing.isNone() ? "no crossing" : "crossing"));
  if (crossing.isNone()) {
    return false;
  }
  exit = crossing;
  return true;
}

NavigationTarget CylinderNavigationPolicy::findExit(
    const GeometryContext& gctx, const VolumeExitArguments& args,
    bool dipped) const {
  using enum CylinderVolumeBounds::Face;

  // Look back along the straight line for the last crossing onto this volume.
  // After a dip, the trajectory came back in through the inner cylinder and
  // left through it right before, within the step.
  const double range =
      dipped ? (args.position - *args.previousPosition).norm() : 0.;

  NavigationTarget last = NavigationTarget::None();
  NavigationTarget secondToLast = NavigationTarget::None();
  for (const Portal* portal : m_portals) {
    if (portal == nullptr ||
        (dipped && portal != m_portals.at(toUnderlying(InnerCylinder)))) {
      continue;
    }
    const MultiIntersection3D intersections = portal->surface().intersect(
        gctx, args.position, args.direction, BoundaryTolerance::None(),
        s_onSurfaceTolerance);
    for (const auto [index, intersection] : enumerate(intersections)) {
      if (intersection.status() != IntersectionStatus::reachable ||
          intersection.pathLength() >= 0. ||
          (dipped && -intersection.pathLength() > range) ||
          !onFace(toLocal(gctx, intersection.position()), args.tolerance)) {
        continue;
      }
      NavigationTarget crossing(intersection,
                                static_cast<IntersectionIndex>(index), *portal,
                                BoundaryTolerance::None());
      if (last.isNone() || crossing.pathLength() > last.pathLength()) {
        secondToLast = last;
        last = crossing;
      } else if (secondToLast.isNone() ||
                 crossing.pathLength() > secondToLast.pathLength()) {
        secondToLast = crossing;
      }
    }
  }

  return dipped ? secondToLast : last;
}
}  // namespace Acts
