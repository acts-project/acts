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

CylinderNavigationPolicy::CylinderNavigationPolicy(const GeometryContext& gctx,
                                                   const TrackingVolume& volume,
                                                   const Logger& logger)
    : m_volume(&volume) {
  ACTS_VERBOSE("CylinderNavigationPolicy constructor for volume "
               << volume.volumeName());
  using enum CylinderVolumeBounds::Face;

  if (m_volume->volumeBounds().type() != VolumeBounds::eCylinder) {
    ACTS_ERROR("CylinderNavigationPolicy can only be used with "
               << "cylinder volumes");
    throw std::invalid_argument(
        "CylinderNavigationPolicy can only be used with "
        "cylinder volumes");
  }

  auto& bounds =
      dynamic_cast<const CylinderVolumeBounds&>(m_volume->volumeBounds());

  double rMin = bounds.get(CylinderVolumeBounds::eMinR);
  m_rMin2 = rMin * rMin;
  double rMax = bounds.get(CylinderVolumeBounds::eMaxR);
  m_rMax2 = rMax * rMax;
  m_halfLengthZ = bounds.get(CylinderVolumeBounds::eHalfLengthZ);

  if (rMin == 0) {
    ACTS_ERROR("CylinderNavigationPolicy can only be used with "
               << "non-zero inner radius, which " << m_volume->volumeName()
               << " has not");
    throw std::invalid_argument(
        "CylinderNavigationPolicy can only be used with "
        "non-zero inner radius");
  }

  if (m_volume->portals().size() != 4) {
    ACTS_ERROR("CylinderNavigationPolicy can only be used with "
               << "volumes with 4 portals");
    throw std::invalid_argument(
        "CylinderNavigationPolicy can only be used with "
        "volumes with 4 portals");
  }

  ACTS_VERBOSE("CylinderNavigationPolicy created for volume "
               << volume.volumeName());

  if (!volume.transform().linear().isApprox(SquareMatrix3::Identity())) {
    m_itransform = volume.transform().inverse();
  }

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
      Transform3 localTransform = m_volume->transform().inverse() *
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

  ACTS_VERBOSE("Portal assignment:");
  for (const auto& [i, portal] : enumerate(m_portals)) {
    auto face = static_cast<CylinderVolumeBounds::Face>(i);

    if (portal == nullptr) {
      ACTS_ERROR("Have no portal for " << faceName(face));
      throw std::invalid_argument("Have no portal for " +
                                  std::string{faceName(face)});
    }

    ACTS_VERBOSE("  " << faceName(face) << " -> "
                      << portal->surface().toStream(gctx));
  }
}

void CylinderNavigationPolicy::initializeCandidates(
    [[maybe_unused]] const GeometryContext& gctx,
    const NavigationArguments& args, AppendOnlyNavigationStream& stream,
    const Logger& logger) const {
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
    pos = args.position - m_volume->transform().translation();
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
}  // namespace Acts
