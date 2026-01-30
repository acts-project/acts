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
  // Ensure that the vector of portals is large enough
  if (portalIdx >= m_portalPlacements.size()) {
    m_portalPlacements.resize(portalIdx + 1);
  }
  // Fetch the registered placement
  std::unique_ptr<detail::PortalPlacement>& placement =
      m_portalPlacements[portalIdx];

  // If there is none just create it
  if (!placement) {
    placement.reset(new detail::PortalPlacement(portalIdx, portalToVolTrf, this,
                                                std::move(portalSurface)));
  }

  // Return the cached shared ptr
  return placement->getSharedPtr();
}

const detail::PortalPlacement* VolumePlacementBase::portalPlacement(
    const std::size_t portalIdx) const {
  return portalIdx < m_portalPlacements.size()
             ? m_portalPlacements[portalIdx].get()
             : nullptr;
}

std::size_t VolumePlacementBase::nPortalPlacements() const {
  return m_portalPlacements.size();
}

void VolumePlacementBase::populateContextWithPortals(
    GeometryContext& gctx) const {
  for (const auto& portalPlacement : m_portalPlacements) {
    cachePortalTransform(gctx, portalPlacement->index(),
                         portalPlacement->assembleFullTransform(gctx));
  }
}

}  // namespace Acts
