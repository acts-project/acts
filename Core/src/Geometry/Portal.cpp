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
      (m_volumeLinks[indexFromDirection(NavigationDirection::Backward)]
           .implementation == nullptr)
          ? NavigationDirection::Forward
          : NavigationDirection::Backward;

  if (m_volumeLinks[indexFromDirection(tDir)].implementation == nullptr) {
    throw std::invalid_argument(
        "Portal: trying to fuse portal (keep) with no links.");
  }
  // And now check other direction
  NavigationDirection oDir = (tDir == NavigationDirection::Forward)
                                 ? NavigationDirection::Backward
                                 : NavigationDirection::Forward;

  if (other->m_volumeLinks[indexFromDirection(oDir)].implementation ==
      nullptr) {
    throw std::runtime_error(
        "Portal: trying to fuse portal (waste) with no links.");
  }

  auto odx = indexFromDirection(oDir);
  m_volumeLinks[odx] = std::move(other->m_volumeLinks[odx]);
  m_attachedVolumes[odx] = other->m_attachedVolumes[odx];

  // And finally overwrite
  other = getSharedPtr();
}

void Acts::Experimental::Portal::updateVolumeLink(
    NavigationDirection nDir, ManagedDetectorVolumeLink&& dVolumeLink,
    const std::vector<std::shared_ptr<DetectorVolume>>& attachedVolumes) {
  auto idx = indexFromDirection(nDir);
  m_volumeLinks[idx] = std::move(dVolumeLink);
  m_attachedVolumes[idx] = attachedVolumes;
}

const Acts::Experimental::DetectorVolume*
Acts::Experimental::Portal::nextVolume(const GeometryContext& gctx,
                                       const Vector3& position,
                                       const Vector3& direction) const {
  const Vector3 normal = surface().normal(gctx, position);
  NavigationDirection nDir = (normal.dot(direction) < 0.)
                                 ? NavigationDirection::Backward
                                 : NavigationDirection::Forward;
  const auto& dLink = m_volumeLinks[indexFromDirection(nDir)].delegate;
  if (dLink.connected()) {
    return dLink(gctx, position, direction);
  }
  return nullptr;
}
