// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/Portal.hpp"

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

void Acts::Experimental::Portal::assignSurfaceMaterial(
    std::shared_ptr<const ISurfaceMaterial> material) {
  return m_surface->assignSurfaceMaterial(material);
}

void Acts::Experimental::Portal::assignGeometryId(
    const GeometryIdentifier& geometryId) {
  m_surface->assignGeometryId(geometryId);
}

void Acts::Experimental::Portal::updateVolumeLink(
    NavigationDirection nDir, const DetectorVolumeLink& dVolumeLink,
    DetectorVolumeLinkStore dVolumeLinkStore, bool creationVolume) {
  if (creationVolume) {
    m_creationVolumeDir = nDir;
  }
  if (dVolumeLinkStore != nullptr) {
    m_linkStore.insert(dVolumeLinkStore);
  }
  if (nDir == NavigationDirection::Forward) {
    m_forwardLink = dVolumeLink;
    return;
  }
  m_backwardLink = dVolumeLink;
}

void Acts::Experimental::Portal::updateOutsideVolumeLink(
    const DetectorVolumeLink& dVolumeLink,
    DetectorVolumeLinkStore dVolumeLinkStore, bool overwrite) {
  if (m_creationVolumeDir.has_value()) {
    if (m_creationVolumeDir.value() == NavigationDirection::Forward and
        (not m_backwardLink.connected() or overwrite)) {
      m_backwardLink = dVolumeLink;
    } else if (not m_forwardLink.connected() or overwrite) {
      m_forwardLink = dVolumeLink;
    }
    if (dVolumeLinkStore != nullptr) {
      m_linkStore.insert(dVolumeLinkStore);
    }
  }
}

const Acts::Experimental::DetectorVolume*
Acts::Experimental::Portal::nextVolume(const GeometryContext& gctx,
                                       const Vector3& position,
                                       const Vector3& direction) const {
  const Vector3 normal = surface().normal(gctx, position);
  if (normal.dot(direction) > 0.) {
    return m_forwardLink(gctx, position, direction);
  }
  return m_backwardLink(gctx, position, direction);
}
