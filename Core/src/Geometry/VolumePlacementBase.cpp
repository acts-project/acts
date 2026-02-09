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

VolumePlacementBase::VolumePlacementBase(VolumePlacementBase&& other) noexcept
    : VolumePlacementBase{} {
  (*this) = std::move(other);
}

VolumePlacementBase& VolumePlacementBase::operator=(
    VolumePlacementBase&& other) noexcept {
  if (&other != this) {
    m_portalPlacements = std::move(other.m_portalPlacements);
    for (const auto& portal : m_portalPlacements) {
      portal->m_parent = this;
    }
  }
  return (*this);
}

void VolumePlacementBase::makePortalsAlignable(
    const GeometryContext& gctx, const PortalVec_t& portalsToAlign) {
  if (portalsToAlign.empty()) {
    throw std::invalid_argument(
        "VolumePlacementBase::makePortalsAlignable() - The portals must not be "
        "empty");
  }

  if (m_portalPlacements.size() > 0lu) {
    throw std::runtime_error(
        "VolumePlacementBase::makePortalsAlignable() - Portals were already "
        "registered before");
  }

  for (std::size_t portalIdx = 0lu; portalIdx < nPortalPlacements();
       ++portalIdx) {
    const std::shared_ptr<RegularSurface>& portalSurface =
        portalsToAlign[portalIdx];
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

    m_portalPlacements.emplace_back(new detail::PortalPlacement(
        portalIdx, portalToVolTrf, this, portalSurface));
  }
  // Before leaving ensure that the client knows it needs to reserve
  // space for N portals in the cache backend
  expandTransformCache(gctx, nPortalPlacements());
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

Transform3 VolumePlacementBase::alignPortal(const GeometryContext& gctx,
                                            const std::size_t portalIdx) const {
  assert(portalIdx < m_portalPlacements.size());
  return m_portalPlacements[portalIdx]->assembleFullTransform(gctx);
}

}  // namespace Acts
