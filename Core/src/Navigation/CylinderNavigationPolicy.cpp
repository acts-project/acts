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

CylinderNavigationPolicy::CylinderNavigationPolicy(const GeometryContext& gctx,
                                                   const TrackingVolume& volume,
                                                   const Logger& logger,
                                                   const Config& config)
    : m_cfg(config), m_volume(&volume) {
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
  double rMax = bounds.get(CylinderVolumeBounds::eMaxR);
  double zHalf = bounds.get(CylinderVolumeBounds::eHalfLengthZ);

  if (m_volume->portals().size() != 4) {
    ACTS_ERROR("CylinderNavigationPolicy can only be used with "
               << "volumes with 4 portals");
    throw std::invalid_argument(
        "CylinderNavigationPolicy can only be used with "
        "volumes with 4 portals");
  }

  ACTS_VERBOSE("CylinderNavigationPolicy created for volume "
               << volume.volumeName());
  m_cylinderAxis = volume.transform().rotation() * Vector3::UnitZ();
  ACTS_VERBOSE("Cylinder axis is " << m_cylinderAxis.transpose());

  // Since the volume does not store the shell assignment, we have to recover
  // this from the raw portals
  for (const auto& portal : m_volume->portals()) {
    if (const auto* cylBounds =
            dynamic_cast<const CylinderBounds*>(&portal.surface().bounds());
        cylBounds != nullptr) {
      std::cout << "rMin: " << rMin << " rMax: " << rMax << std::endl;
      std::cout << "cylBounds->get(CylinderBounds::eR): "
                << cylBounds->get(CylinderBounds::eR) << std::endl;
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
      if (std::abs(localZ - zHalf) < s_onSurfaceTolerance) {
        m_portals.at(toUnderlying(PositiveDisc)) = &portal;
        continue;
      }

      if (std::abs(localZ + zHalf) < s_onSurfaceTolerance) {
        m_portals.at(toUnderlying(NegativeDisc)) = &portal;
        continue;
      }
    }
  }

  ACTS_VERBOSE("Portal assignment:");
  for (const auto& [i, portal] : enumerate(m_portals)) {
    auto face = static_cast<CylinderVolumeBounds::Face>(i);

    std::string faceName;
    switch (face) {
      case PositiveDisc:
        faceName = "PositiveDisc";
        break;
      case NegativeDisc:
        faceName = "NegativeDisc";
        break;
      case OuterCylinder:
        faceName = "OuterCylinder";
        break;
      case InnerCylinder:
        faceName = "InnerCylinder";
        break;
      default:
        faceName = "Unknown";
    }
    if (portal == nullptr) {
      ACTS_ERROR("Have no portal for " << faceName);
      throw std::invalid_argument("Have no portal for " + faceName);
    }

    ACTS_VERBOSE("  " << faceName << " -> "
                      << portal->surface().toStream(gctx));
  }
}

void CylinderNavigationPolicy::initializeCandidates(
    const NavigationArguments& args, AppendOnlyNavigationStream& stream,
    const Logger& logger) const {
  Vector3 pos = args.position - m_volume->center();
  Vector3 dir = args.direction;

  double dot = dir.dot(m_cylinderAxis);
}

void CylinderNavigationPolicy::connect(NavigationDelegate& delegate) const {
  connectDefault<CylinderNavigationPolicy>(delegate);
}
}  // namespace Acts
