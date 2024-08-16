// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/PortalShell.hpp"

#include "Acts/Geometry/BoundarySurfaceFace.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts {

SingleCylinderPortalShell::SingleCylinderPortalShell(TrackingVolume& volume) {
  assert(volume.volumeBounds().type() == VolumeBounds::BoundsType::eCylinder);
  const auto& bounds =
      dynamic_cast<const CylinderVolumeBounds&>(volume.volumeBounds());

  std::vector<OrientedSurface> orientedSurfaces =
      bounds.orientedSurfaces(volume.transform());

  auto handle = [&](Face face, std::size_t from) {
    const auto& source = orientedSurfaces.at(from);
    m_portals.at(toUnderlying(face)) =
        std::make_shared<Portal>(source.direction, source.surface, volume);
  };

  if (orientedSurfaces.size() == 6) {
    // Fully equipped cylinder
    handle(PositiveDisc, positiveFaceXY);
    handle(NegativeDisc, negativeFaceXY);
    handle(OuterCylinder, tubeOuterCover);
    handle(InnerCylinder, tubeInnerCover);
    handle(NegativePhiPlane, tubeSectorNegativePhi);
    handle(PositivePhiPlane, tubeSectorPositivePhi);
  } else if (orientedSurfaces.size() == 5) {
    // Phi sector but no inner cylinder (rMin == 0)
    handle(PositiveDisc, positiveFaceXY);
    handle(NegativeDisc, negativeFaceXY);
    handle(OuterCylinder, tubeOuterCover);
    // Skip inner tube cover, requires offsetting
    handle(NegativePhiPlane, tubeSectorNegativePhi - 1);
    handle(PositivePhiPlane, tubeSectorPositivePhi - 1);
  } else if (orientedSurfaces.size() == 4) {
    // No phi sector but rMin > 0
    handle(PositiveDisc, positiveFaceXY);
    handle(NegativeDisc, negativeFaceXY);
    handle(OuterCylinder, tubeOuterCover);
    handle(InnerCylinder, tubeInnerCover);
  } else if (orientedSurfaces.size() == 3) {
    // No phi sector and rMin == 0
    handle(PositiveDisc, positiveFaceXY);
    handle(NegativeDisc, negativeFaceXY);
    handle(OuterCylinder, tubeOuterCover);
  } else {
    throw std::invalid_argument("Invalid number of oriented surfaces");
  }
}

Portal* SingleCylinderPortalShell::portal(Face face) {
  return m_portals.at(toUnderlying(face)).get();
}

std::size_t SingleCylinderPortalShell::size() const {
  std::size_t count = 0;
  std::ranges::for_each(
      m_portals, [&count](const auto& portal) { count += portal ? 1 : 0; });
  return count;
}

}  // namespace Acts
