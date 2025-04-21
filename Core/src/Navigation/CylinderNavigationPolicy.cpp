// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/CylinderNavigationPolicy.hpp"

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
                                                   const Logger& logger,
                                                   const Config& config)
    : m_cfg(config), m_volume(&volume), m_portals{} {
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
               << "non-zero inner radius");
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

  m_itransform = volume.transform().inverse();

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
      Transform3 localTransform =
          m_volume->transform().inverse() * portal.surface().transform(gctx);
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
    const NavigationArguments& args, AppendOnlyNavigationStream& stream,
    const Logger& logger) const {
  using enum CylinderVolumeBounds::Face;
  ACTS_VERBOSE("CylinderNavigationPolicy::initializeCandidates: "
               << "gpos: " << args.position.transpose()
               << " gdir: " << args.direction.transpose());

  const Vector3 pos = m_itransform * args.position;
  const Vector3 dir = m_itransform.linear() * args.direction;

  ACTS_VERBOSE("-> lpos: " << pos.transpose() << " ldir: " << dir.transpose());

  auto add = [&](auto face) {
    ACTS_VERBOSE("~~> Adding portal candidate " << faceName(face));
    stream.addPortalCandidate(*m_portals.at(toUnderlying(face)));
  };

  double rpos2 = pos[0] * pos[0] + pos[1] * pos[1];

  bool hitInner = false;
  // Only run inner cylinder check if we're not ON the inner cylinder
  if (std::abs(rpos2 - m_rMin2) > s_onSurfaceTolerance) {
    // Calculate if we could hit the inner cylinder
    Vector2 dir2 = dir.head<2>().normalized();
    double d = -1 * pos.head<2>().dot(dir2);
    if (d > 0) {  // only check distance if the point of closest approach is in
                  // front of the point
      Vector2 poc = pos.head<2>() + d * dir2;
      double r2 = poc.dot(poc);
      hitInner = r2 < m_rMin2;
      if (hitInner) {
        // Could hit the inner cylinder
        add(InnerCylinder);
      }
    }
  }

  double zDisk = (dir[2] > 0) ? m_halfLengthZ : -m_halfLengthZ;
  if (std::abs(dir[2]) > s_onSurfaceTolerance) {
    // Not parallel to the disc, see if we're inside the disc

    double t = (zDisk - pos[2]) / dir[2];

    Vector3 ix = pos + t * dir;
    assert(std::abs(ix[0] - zDisk) < s_onSurfaceTolerance);
    double r2 = ix[0] * ix[0] + ix[1] * ix[1];
    if (r2 < m_rMax2 && r2 > m_rMin2) {
      // Will hit this disk!
      add(dir[2] > 0 ? PositiveDisc : NegativeDisc);
      // If we hit the disk inside the radius, we can't hit the outer cylinder!
      return;
    }
  }

  if (!hitInner) {
    // If we don't hit the inner cylinder, we can hit the outer cylinder
    add(OuterCylinder);
  }
}

void CylinderNavigationPolicy::connect(NavigationDelegate& delegate) const {
  // @TODO: Implement optimization for shift only transform, where we can skip the rotation
  connectDefault<CylinderNavigationPolicy>(delegate);
}
}  // namespace Acts
