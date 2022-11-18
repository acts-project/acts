// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Portal.hpp"

#include "Acts/Surfaces/Surface.hpp"

Acts::Experimental::Portal::Portal(std::shared_ptr<Surface> surface)
    : m_surface(std::move(surface)) {}

std::shared_ptr<Acts::Experimental::Portal>
Acts::Experimental::Portal::getSharedPtr() {
  return shared_from_this();
}

std::shared_ptr<const Acts::Experimental::Portal>
Acts::Experimental::Portal::getSharedPtr() const {
  return shared_from_this();
}

void Acts::Experimental::Portal::assignGeometryId(
    const GeometryIdentifier& geometryId) {
  m_surface->assignGeometryId(geometryId);
}

void Acts::Experimental::Portal::fuse(std::shared_ptr<Portal>& other) {
  // Determine this directioon
  NavigationDirection tDir =
      (not m_volumeUpdators[indexFromDirection(NavigationDirection::Backward)]
               .connected())
          ? NavigationDirection::Forward
          : NavigationDirection::Backward;

  if (not m_volumeUpdators[indexFromDirection(tDir)].connected()) {
    throw std::invalid_argument(
        "Portal: trying to fuse portal (keep) with no links.");
  }
  // And now check other direction
  NavigationDirection oDir = invertDirection(tDir);
  if (not other->m_volumeUpdators[indexFromDirection(oDir)].connected()) {
    throw std::runtime_error(
        "Portal: trying to fuse portal (waste) with no links.");
  }

  auto odx = indexFromDirection(oDir);
  m_volumeUpdators[odx] = std::move(other->m_volumeUpdators[odx]);
  m_attachedVolumes[odx] = other->m_attachedVolumes[odx];
  // And finally overwrite the original portal
  other = getSharedPtr();
}

void Acts::Experimental::Portal::assignDetectorVolumeUpdator(
    NavigationDirection nDir, DetectorVolumeUpdator&& dVolumeUpdator,
    const std::vector<std::shared_ptr<DetectorVolume>>& attachedVolumes) {
  auto idx = indexFromDirection(nDir);
  m_volumeUpdators[idx] = std::move(dVolumeUpdator);
  m_attachedVolumes[idx] = attachedVolumes;
}

void Acts::Experimental::Portal::assignDetectorVolumeUpdator(
    DetectorVolumeUpdator&& dVolumeUpdator,
    const std::vector<std::shared_ptr<DetectorVolume>>& attachedVolumes) {
  // Check and throw exceptions
  if (not m_volumeUpdators[0u].connected() and
      not m_volumeUpdators[1u].connected()) {
    throw std::runtime_error("Portal: portal has no link on either side.");
  }
  if (m_volumeUpdators[0u].connected() and m_volumeUpdators[1u].connected()) {
    throw std::runtime_error("Portal: portal already has links on both sides.");
  }
  size_t idx = m_volumeUpdators[0u].connected() ? 1u : 0u;
  m_volumeUpdators[idx] = std::move(dVolumeUpdator);
  m_attachedVolumes[idx] = attachedVolumes;
}

void Acts::Experimental::Portal::updateDetectorVolume(
    const GeometryContext& gctx, NavigationState& nState) const {
  const auto& position = nState.position;
  const auto& direction = nState.direction;
  const Vector3 normal = surface().normal(gctx, position);
  NavigationDirection nDir = directionFromStepSize(normal.dot(direction));
  const auto& vUpdator = m_volumeUpdators[indexFromDirection(nDir)];
  if (vUpdator.connected()) {
    vUpdator(gctx, nState);
  } else {
    nState.currentVolume = nullptr;
  }
}
