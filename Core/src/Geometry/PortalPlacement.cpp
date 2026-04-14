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

PortalPlacement::PortalPlacement(ConstructTag /*tag*/,
                                 const std::size_t portalIdx,
                                 const Transform3& portalTrf,
                                 const VolumePlacementBase* parent)
    : m_portalToVolumeCenter{portalTrf},
      m_parent{parent},
      m_portalIdx{portalIdx} {}

Transform3 PortalPlacement::assembleFullTransform(
    const GeometryContext& gctx) const {
  return m_parent->localToGlobalTransform(gctx) * portalToVolumeCenter();
}

const Transform3& PortalPlacement::localToGlobalTransform(
    const GeometryContext& gctx) const {
  return m_parent->portalLocalToGlobal(gctx, m_portalIdx);
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

}  // namespace Acts::detail
