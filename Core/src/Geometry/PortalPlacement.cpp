// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/detail/PortalPlacement.hpp"

#include "Acts/Geometry/VolumePlacementBase.hpp"

namespace Acts::detail {

PortalPlacement::PortalPlacement(const std::size_t portalIdx,
                                 const Transform3& portalTrf,
                                 VolumePlacementBase* parent,
                                 std::shared_ptr<RegularSurface> surface)
    : m_portalToVolumeCenter{portalTrf},
      m_surface{std::move(surface)},
      m_parent{parent},
      m_portalIdx{portalIdx} {
  assert(m_surface != nullptr);
  m_surface->assignSurfacePlacement(*this);
}

Transform3 PortalPlacement::assembleFullTransform(
    const GeometryContext& gctx) const {
  return m_parent->localToGlobalTransform(gctx) * portalToVolumeCenter();
}

const Transform3& PortalPlacement::localToGlobalTransform(
    const GeometryContext& gctx) const {
  return m_parent->portalLocalToGlobal(gctx, m_portalIdx);
}

const RegularSurface& PortalPlacement::surface() const {
  return *m_surface;
}

RegularSurface& PortalPlacement::surface() {
  return *m_surface;
}

bool PortalPlacement::isSensitive() const {
  return false;
}

std::size_t PortalPlacement::index() const {
  return m_portalIdx;
}

const Transform3& PortalPlacement::portalToVolumeCenter() const {
  return m_portalToVolumeCenter;
}

std::shared_ptr<RegularSurface> PortalPlacement::surfacePtr() {
  return m_surface;
}

std::shared_ptr<const RegularSurface> PortalPlacement::surfacePtr() const {
  return m_surface;
}
}  // namespace Acts::detail
