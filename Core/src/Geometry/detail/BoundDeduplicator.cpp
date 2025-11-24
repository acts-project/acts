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

#define IMPL_SURF_DEDUPLICATION(ENUM_CASE, SURFACE_T)   \
  case ENUM_CASE: {                                     \
    auto& castSurf = dynamic_cast<SURFACE_T&>(surface); \
    if (castSurf.boundsPtr()) {                         \
      castSurf.assignSurfaceBounds(                     \
          m_surfFactory.insert(castSurf.boundsPtr()));  \
    }                                                   \
    break;                                              \
  }

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
    IMPL_SURF_DEDUPLICATION(Cone, ConeSurface);
    IMPL_SURF_DEDUPLICATION(Cylinder, CylinderSurface);
    IMPL_SURF_DEDUPLICATION(Disc, DiscSurface);
    IMPL_SURF_DEDUPLICATION(Plane, PlaneSurface);
    IMPL_SURF_DEDUPLICATION(Straw, LineSurface);
    default:
      throw std::invalid_argument(std::format(
          "BoundDeduplicator::visitSurface() - The surface {:} is not yet "
          "supported",
          Surface::s_surfaceTypeNames[toUnderlying(surface.type())]));
  }
}

void BoundDeduplicator::visitBoundarySurface(
    BoundarySurfaceT<TrackingVolume>& boundary) {
  visitSurface(boundary.surfaceRepresentation());
}

}  // namespace Acts::detail
