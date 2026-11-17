// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/detail/CylinderPortalBounds.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace Acts::detail {

CylinderPortalBounds::CylinderPortalBounds(const GeometryContext& gctx,
                                           const TrackingVolume& volume) {
  if (volume.volumeBounds().type() != VolumeBounds::eCylinder ||
      (volume.portals().size() != 4 && volume.portals().size() != 6)) {
    return;
  }
  const auto& bounds =
      static_cast<const CylinderVolumeBounds&>(volume.volumeBounds());
  nestedContainer = volume.portals().size() == 6;
  if ((!nestedContainer && bounds.get(CylinderVolumeBounds::eMinR) == 0.) ||
      bounds.get(CylinderVolumeBounds::eHalfPhiSector) < std::numbers::pi) {
    return;
  }
  const auto isAxial = [](const Transform3& transform) {
    const Vector3 axis = transform.linear().col(2);
    return axis.head<2>().squaredNorm() <= 1e-24 &&
           std::abs(std::abs(axis.z()) - 1.) <= 1e-12;
  };
  if (!isAxial(volume.localToGlobalTransform(gctx))) {
    return;
  }
  center = volume.center(gctx).head<2>();
  // Validate every supporting face before publishing any descriptor. Portal
  // count alone does not establish the topology, and shared surfaces need not
  // have the owner's transform. These strict tolerances only identify coaxial
  // geometry; they are not tracking or boundary tolerances.
  std::size_t cylinderCount = 0;
  std::size_t discCount = 0;
  for (const auto& portal : volume.portals()) {
    const auto& surface = portal.surface();
    if (!isAxial(surface.localToGlobalTransform(gctx)) ||
        (surface.center(gctx).head<2>() - center).squaredNorm() > 1e-24) {
      return;
    }
    if (surface.type() == Surface::Cylinder) {
      const auto* cylinder =
          dynamic_cast<const CylinderBounds*>(&surface.bounds());
      if (cylinder == nullptr ||
          cylinder->get(CylinderBounds::eHalfPhiSector) < std::numbers::pi) {
        return;
      }
      ++cylinderCount;
    } else if (surface.type() == Surface::Disc) {
      const auto* disc = dynamic_cast<const DiscBounds*>(&surface.bounds());
      if (disc == nullptr || !disc->coversFullAzimuth()) {
        return;
      }
      ++discCount;
    } else {
      return;
    }
  }
  if (cylinderCount != 2 || discCount != (nestedContainer ? 4u : 2u)) {
    return;
  }
  lowerZ = std::numeric_limits<double>::infinity();
  upperZ = -std::numeric_limits<double>::infinity();
  double innerRadius = bounds.get(CylinderVolumeBounds::eMinR);
  if (nestedContainer) {
    // A solid parent with one nested cylindrical child has two cylindrical
    // interfaces. While in its radial annulus, the child outer interface
    // acts as the parent's inward exit. Preserve the native Portal resolver.
    // This deliberately covers only the exported coaxial two-cylinder case.
    innerRadius = std::numeric_limits<double>::infinity();
    for (const auto& portal : volume.portals()) {
      const auto& surface = portal.surface();
      if (surface.type() != Surface::Cylinder) {
        continue;
      }
      const auto& cylinder =
          static_cast<const CylinderBounds&>(surface.bounds());
      innerRadius = std::min(innerRadius, cylinder.get(CylinderBounds::eR));
    }
    if (!(innerRadius > 0.) ||
        !(innerRadius < bounds.get(CylinderVolumeBounds::eMaxR))) {
      return;
    }
  }
  const double outerRadius = bounds.get(CylinderVolumeBounds::eMaxR);
  innerRadius2 = innerRadius * innerRadius;
  outerRadius2 = outerRadius * outerRadius;
  for (const auto& portal : volume.portals()) {
    const auto& surface = portal.surface();
    if (surface.type() == Surface::Cylinder) {
      const auto& cylinder =
          static_cast<const CylinderBounds&>(surface.bounds());
      const double r = cylinder.get(CylinderBounds::eR);
      if (std::abs(r - innerRadius) < s_onSurfaceTolerance) {
        innerCylinder = &portal;
        innerMinZ = surface.center(gctx).z() -
                    cylinder.get(CylinderBounds::eHalfLengthZ);
        innerMaxZ = surface.center(gctx).z() +
                    cylinder.get(CylinderBounds::eHalfLengthZ);
      } else if (std::abs(r - outerRadius) < s_onSurfaceTolerance) {
        outerCylinder = &portal;
      }
      continue;
    }
    const auto& transform = surface.localToGlobalTransform(gctx);
    const double z = transform.translation().z();
    if (nestedContainer && containerDiscCount < containerDiscs.size()) {
      containerDiscs[containerDiscCount++] = {&portal, z};
    }
    if (z < lowerZ) {
      lowerZ = z;
      lowerDisc = &portal;
    }
    if (z > upperZ) {
      upperZ = z;
      upperDisc = &portal;
    }
  }
  cachedVolume = &volume;
  if (lowerDisc == upperDisc) {
    lowerDisc = nullptr;
    upperDisc = nullptr;
  }
}

bool CylinderPortalBounds::valid() const {
  return cachedVolume != nullptr && innerCylinder != nullptr &&
         outerCylinder != nullptr && lowerDisc != nullptr &&
         upperDisc != nullptr;
}

void CylinderPortalBounds::attach(AppendOnlyNavigationStream& stream) const {
  if (valid()) {
    stream.m_stream->m_cylinderPortalBounds = this;
  }
}

}  // namespace Acts::detail
