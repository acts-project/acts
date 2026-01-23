// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/VolumePlacementBase.hpp"

#include <cassert>
namespace Acts {

std::shared_ptr<RegularSurface> VolumePlacementBase::makePortalAlignable(
    const std::size_t portalIdx, const Transform3& portalToVolTrf,
    std::shared_ptr<RegularSurface>&& portalSurface) {
  if (portalIdx >= m_portalPlacements.size()) {
    m_portalPlacements.resize(portalIdx + 1);
  }
  std::unique_ptr<detail::PortalPlacement>& placement =
      m_portalPlacements[portalIdx];
  if (!placement) {
    placement.reset(new detail::PortalPlacement(portalIdx, portalToVolTrf, this,
                                                std::move(portalSurface)));
  }
  assert(portalToVolTrf.isApprox(placement->portalToVolumeCenter()));
  return dynamic_pointer_cast<RegularSurface>(
      placement->surface().getSharedPtr());
}

}  // namespace Acts
