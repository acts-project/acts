// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Portal.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>

namespace Acts::Experimental {

Portal::Portal(std::shared_ptr<RegularSurface> surface)
    : m_surface(std::move(surface)) {
  throw_assert(m_surface, "Portal surface is nullptr");
}

const Acts::RegularSurface& Portal::surface() const {
  return *m_surface.get();
}

Acts::RegularSurface& Portal::surface() {
  return *m_surface.get();
}

const Portal::DetectorVolumeUpdaters& Portal::detectorVolumeUpdaters() const {
  return m_volumeUpdaters;
}

Portal::AttachedDetectorVolumes& Portal::attachedDetectorVolumes() {
  return m_attachedVolumes;
}

void Portal::assignGeometryId(const GeometryIdentifier& geometryId) {
  m_surface->assignGeometryId(geometryId);
}

std::shared_ptr<Portal> Portal::fuse(std::shared_ptr<Portal>& aPortal,
                                     std::shared_ptr<Portal>& bPortal) {
  if (aPortal == bPortal) {
    return aPortal;
  }

  auto bothConnected = [](const auto& p) {
    return p.m_volumeUpdaters[0u].connected() &&
           p.m_volumeUpdaters[1u].connected();
  };

  auto noneConnected = [](const auto& p) {
    return !p.m_volumeUpdaters[0u].connected() &&
           !p.m_volumeUpdaters[1u].connected();
  };

  if (bothConnected(*aPortal) || bothConnected(*bPortal)) {
    throw std::invalid_argument(
        "Portal: trying to fuse two portals where at least one has links on "
        "both sides.");
  }

  if (noneConnected(*aPortal) || noneConnected(*bPortal)) {
    throw std::invalid_argument(
        "Portal: trying to fuse two portals where at least on has no links.");
  }

  // We checked they're not both empty, so one of them must be connected
  Direction aDir = (aPortal->m_volumeUpdaters[0].connected())
                       ? Direction::fromIndex(0)
                       : Direction::fromIndex(1);
  Direction bDir = aDir.invert();

  // And now check other direction
  if (!bPortal->m_volumeUpdaters[bDir.index()].connected()) {
    throw std::runtime_error(
        "Portal: trying to fuse portal (discard) with no links.");
  }

  const auto& aSurface = aPortal->surface();
  const auto& bSurface = bPortal->surface();

  // @TODO: There's no safety against fusing portals with different surfaces
  std::shared_ptr<Portal> fused = std::make_shared<Portal>(aPortal->m_surface);

  if (aSurface.surfaceMaterial() != nullptr &&
      bSurface.surfaceMaterial() != nullptr) {
    throw std::runtime_error(
        "Portal: both surfaces have surface material, fusing will lead to "
        "information loss.");
  } else if (aSurface.surfaceMaterial() != nullptr) {
    fused->m_surface = aPortal->m_surface;
  } else if (bSurface.surfaceMaterial() != nullptr) {
    fused->m_surface = bPortal->m_surface;
  }

  fused->m_volumeUpdaters[aDir.index()] =
      std::move(aPortal->m_volumeUpdaters[aDir.index()]);
  fused->m_attachedVolumes[aDir.index()] =
      std::move(aPortal->m_attachedVolumes[aDir.index()]);

  fused->m_volumeUpdaters[bDir.index()] =
      std::move(bPortal->m_volumeUpdaters[bDir.index()]);
  fused->m_attachedVolumes[bDir.index()] =
      std::move(bPortal->m_attachedVolumes[bDir.index()]);

  return fused;
}

void Portal::assignDetectorVolumeUpdater(
    Direction dir, DetectorVolumeUpdater dVolumeUpdater,
    std::vector<std::shared_ptr<DetectorVolume>> attachedVolumes) {
  auto idx = dir.index();
  m_volumeUpdaters[idx] = std::move(dVolumeUpdater);
  m_attachedVolumes[idx] = std::move(attachedVolumes);
}

void Portal::assignDetectorVolumeUpdater(
    DetectorVolumeUpdater dVolumeUpdater,
    std::vector<std::shared_ptr<DetectorVolume>> attachedVolumes) {
  // Check and throw exceptions
  if (!m_volumeUpdaters[0u].connected() && !m_volumeUpdaters[1u].connected()) {
    throw std::runtime_error("Portal: portal has no link on either side.");
  }
  if (m_volumeUpdaters[0u].connected() && m_volumeUpdaters[1u].connected()) {
    throw std::runtime_error("Portal: portal already has links on both sides.");
  }
  std::size_t idx = m_volumeUpdaters[0u].connected() ? 1u : 0u;
  m_volumeUpdaters[idx] = std::move(dVolumeUpdater);
  m_attachedVolumes[idx] = std::move(attachedVolumes);
}

void Portal::updateDetectorVolume(const GeometryContext& gctx,
                                  NavigationState& nState) const {
  const auto& position = nState.position;
  const auto& direction = nState.direction;
  const Vector3 normal = surface().normal(gctx, position);
  Direction dir = Direction::fromScalar(normal.dot(direction));
  const auto& vUpdater = m_volumeUpdaters[dir.index()];
  if (vUpdater.connected()) {
    vUpdater(gctx, nState);
  } else {
    nState.currentVolume = nullptr;
  }
}

}  // namespace Acts::Experimental
