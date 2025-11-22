// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/detail/BoundDeduplicator.hpp"

#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"

#include <format>

namespace Acts::detail {
void BoundDeduplicator::visitVolume(TrackingVolume& volume) {
  volume.assignVolumeBounds(m_volFactory.insert(volume.volumeBoundsPtr()));
}
void BoundDeduplicator::visitPortal(Portal& portal) {
  visitSurface(portal.surface());
}

void BoundDeduplicator::visitSurface(Surface& surface) {
  switch (surface.type()) {
    using enum Surface::SurfaceType;
    case Cone: {
      auto& castSurf = static_cast<ConeSurface&>(surface);
      castSurf.m_bounds = m_surfFactory.insert(castSurf.m_bounds);
      break;
    }
    case Cylinder: {
      auto& castSurf = static_cast<CylinderSurface&>(surface);
      castSurf.m_bounds = m_surfFactory.insert(castSurf.m_bounds);
      break;
    }
    case Disc: {
      auto& castSurf = static_cast<CylinderSurface&>(surface);
      castSurf.m_bounds = m_surfFactory.insert(castSurf.m_bounds);
      break;
    }
    case Plane: {
      auto& castSurf = static_cast<PlaneSurface&>(surface);
      castSurf.m_bounds = m_surfFactory.insert(castSurf.m_bounds);
      break;
    }
    case Straw: {
      auto& castSurf = static_cast<LineSurface&>(surface);
      castSurf.m_bounds = m_surfFactory.insert(castSurf.m_bounds);
      break;
    }
    default:
      throw std::invalid_argument(std::format(
          "BoundDeduplicator::visitSurface() - The surface {:} is not yet "
          "supported",
          Surface::s_surfaceTypeNames[toUnderlying(surface.type())]));
  }
}

void BoundDeduplicator::visitBoundarySurface(
    BoundarySurfaceT<TrackingVolume>& boundary) {
  visitSurface(const_cast<RegularSurface&>(boundary.surfaceRepresentation()));
}
}  // namespace Acts::detail
