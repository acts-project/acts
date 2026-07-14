// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/detail/PortalShellFactory.hpp"

#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CylinderPortalShell.hpp"
#include "Acts/Geometry/DiamondPortalShell.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidPortalShell.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <stdexcept>

namespace Acts::detail {

std::unique_ptr<PortalShellBase> makeSinglePortalShell(
    const GeometryContext& gctx, TrackingVolume& volume) {
  switch (volume.volumeBounds().type()) {
    case VolumeBounds::eCylinder:
      return std::make_unique<SingleCylinderPortalShell>(gctx, volume);
    case VolumeBounds::eCuboid:
      return std::make_unique<SingleCuboidPortalShell>(gctx, volume);
    case VolumeBounds::eTrapezoid:
      return std::make_unique<SingleTrapezoidPortalShell>(gctx, volume);
    case VolumeBounds::eDiamond:
      return std::make_unique<SingleDiamondPortalShell>(gctx, volume);
    default:
      throw std::logic_error("Volume type is not supported");
  }
}

}  // namespace Acts::detail
