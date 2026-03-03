// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/VolumePlacementBase.hpp"

#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <format>
#include <sstream>
namespace Acts {

VolumePlacementBase::VolumePlacementBase() noexcept = default;

VolumePlacementBase::~VolumePlacementBase() = default;

void VolumePlacementBase::makePortalsAlignable(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<RegularSurface>>& portalsToAlign) {
  if (portalsToAlign.empty()) {
    throw std::invalid_argument(
        "VolumePlacementBase::makePortalsAlignable() - The portals must not be "
        "empty");
  }

  if (!m_portalPlacements.empty()) {
    throw std::runtime_error(
        "VolumePlacementBase::makePortalsAlignable() - Portals were already "
        "registered before");
  }

  for (const auto& [portalIdx, portalSurface] : enumerate(portalsToAlign)) {
    if (portalSurface->surfacePlacement() != nullptr) {
      throw std::invalid_argument(std::format(
          "VolumePlacementBase::makePortalsAlignable() - The {:}-th surface is "
          "already connected to the alignment system",
          portalIdx));
    }
    // Calculate the portal transform w.r.t the current alignment
    const Transform3 portalToVolTrf =
        globalToLocalTransform(gctx) *
        portalSurface->localToGlobalTransform(gctx);

    m_portalPlacements.emplace_back(std::make_unique<detail::PortalPlacement>(
        portalIdx, portalToVolTrf, this, portalSurface));
  }
}

const detail::PortalPlacement* VolumePlacementBase::portalPlacement(
    const std::size_t portalIdx) const {
  return m_portalPlacements.at(portalIdx).get();
}

detail::PortalPlacement* VolumePlacementBase::portalPlacement(
    const std::size_t portalIdx) {
  return m_portalPlacements.at(portalIdx).get();
}

std::size_t VolumePlacementBase::nPortalPlacements() const {
  return m_portalPlacements.size();
}

Transform3 VolumePlacementBase::alignPortal(const GeometryContext& gctx,
                                            const std::size_t portalIdx) const {
  return m_portalPlacements.at(portalIdx)->assembleFullTransform(gctx);
}

}  // namespace Acts
