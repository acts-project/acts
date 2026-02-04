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

VolumePlacementBase::PortalVec_t VolumePlacementBase::makePortalsAlignable(
    PortalVec_t portalsToAlign) {
  // It's safe to use the default geometry context as all incoming portal
  // surfaces are explicitly requested to have no surface placement and hence
  // the cached transform is returned
  const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  const bool callExpand = portalsToAlign.size() != m_portalPlacements.size();

  m_portalPlacements.resize(
      std::max(portalsToAlign.size(), m_portalPlacements.size()));
  if (m_portalPlacements.size() != portalsToAlign.size()) {
    throw std::logic_error(
        std::format("VolumePlacementBase::makePortalsAlignable() - Has been "
                    "called previously for a different type of volume bounds. "
                    "Already registered portals {:} vs incoming portals {:}",
                    m_portalPlacements.size(), portalsToAlign.size()));
  }
  for (std::size_t portalIdx = 0lu; portalIdx < nPortalPlacements();
       ++portalIdx) {
    std::shared_ptr<RegularSurface>& portalSurface =
        portalsToAlign[portalIdx].surface;
    if (portalSurface->surfacePlacement() != nullptr) {
      throw std::invalid_argument(std::format(
          "VolumePlacementBase::makePortalsAlignable() - The {:}-th surface is "
          "already connected to the alignment system",
          portalIdx));
    }
    // It is expected that the volume bounds are calling the orientedSurfaces
    // method with an identity transform
    const Transform3 portalToVolTrf =
        portalSurface->localToGlobalTransform(gctx);
    std::unique_ptr<detail::PortalPlacement>& placement =
        m_portalPlacements[portalIdx];
    // Just create a new placement which makes the surface alignable
    if (!placement) {
      placement.reset(new detail::PortalPlacement(portalIdx, portalToVolTrf,
                                                  this, portalSurface));
    } else {
      // The method has been called previously and it needs to be checked that
      // the i-th surface is at the same place as the already registered surface
      if (!portalToVolTrf.isApprox(placement->portalToVolumeCenter())) {
        throw std::invalid_argument(
            std::format("VolumePlacementBase::makePortalsAlignable() -  The "
                        "{:}-th surface does not align with the one from a  "
                        "previous call:\n -- exist {:}\n -- this call {:}",
                        portalIdx, toString(placement->portalToVolumeCenter()),
                        toString(portalToVolTrf)));
      }
      if (!(placement->surface().bounds() == portalSurface->bounds())) {
        std::stringstream sstr{};
        sstr << "VolumePlacementBase::makePortalsAlignable() -  The bounds for "
                "the ";
        sstr << portalIdx
             << "-th surface differ with the ones from a previous call: \n";
        sstr << " -- exist: " << placement->surface().bounds() << "\n";
        sstr << " -- this call: " << portalSurface->bounds();
        throw std::invalid_argument(sstr.str());
      }
      // overwrite the pointer
      portalSurface = placement->surfacePtr();
    }
  }
  // Before leaving ensure that the client knows it needs to reserve
  // space for N portals in the cache backend
  if (callExpand) {
    expandTransformCache(portalsToAlign.size());
  }
  return portalsToAlign;
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
