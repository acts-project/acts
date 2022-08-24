// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/DetectorVolume.hpp"

#include "Acts/Experimental/NavigationState.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"

Acts::Experimental::DetectorVolume::DetectorVolume(
    const GeometryContext& gctx, const Transform3& transform,
    std::unique_ptr<VolumeBounds> bounds,
    const PortalGenerator& portalGenerator, const std::string& name)
    : m_transform(transform),
      m_bounds(std::move(bounds)),
      m_volumeMaterial(nullptr),
      m_name(name) {
  if (m_bounds == nullptr) {
    throw std::invalid_argument(
        "DetectorVolume: construction with nullptr bounds");
  }

  // Contruct the portals
  construct(gctx, portalGenerator);
}

void Acts::Experimental::DetectorVolume::construct(
    const GeometryContext& gctx, const PortalGenerator& portalGenerator) {
  // Create portals with the given generator
  auto portalSurfaces =
      portalGenerator(transform(gctx), *(m_bounds.get()), *this);
  m_portals = ObjectStore<std::shared_ptr<Portal>>(portalSurfaces);
}

std::shared_ptr<Acts::Experimental::DetectorVolume>
Acts::Experimental::DetectorVolume::getSharedPtr() {
  return shared_from_this();
}

std::shared_ptr<const Acts::Experimental::DetectorVolume>
Acts::Experimental::DetectorVolume::getSharedPtr() const {
  return shared_from_this();
}

bool Acts::Experimental::DetectorVolume::inside(const GeometryContext& gctx,
                                                const Vector3& position,
                                                ActsScalar tolerance) const {
  Vector3 posInVolFrame((transform(gctx).inverse()) * position);
  return (volumeBounds()).inside(posInVolFrame, tolerance);
}

void Acts::Experimental::DetectorVolume::updateNavigationStatus(
    NavigationState& nState, const GeometryContext& gctx,
    const Vector3& position, const Vector3& direction, ActsScalar absMomentum,
    ActsScalar charge) const {
  m_navigationStateUpdator(nState, gctx, position, direction, absMomentum,
                           charge);
  return;
}

void Acts::Experimental::DetectorVolume::resize(
    const GeometryContext& gctx, std::unique_ptr<VolumeBounds> rBounds,
    const PortalGenerator& portalGenerator) {
  if (rBounds == nullptr or rBounds->type() != m_bounds->type()) {
    throw std::invalid_argument(
        "DetectorVolume: wrong bound type provided for resize(..) call");
  }
  m_bounds = std::move(rBounds);
  construct(gctx, portalGenerator);
}

void Acts::Experimental::DetectorVolume::lock(
    const GeometryIdentifier& geometryId) {
  m_geometryId = geometryId;

  // Assign the boundary Identifier
  GeometryIdentifier portalId = geometryId;
  for (auto [i, p] : enumerate(m_portals.internal)) {
    portalId.setBoundary(i + 1);
    p->assignGeometryId(portalId);
  }

  // Assign the sensitive/surface Identifier
  GeometryIdentifier sensitiveId = geometryId;
  /// @todo add passive count
  for (auto [i, s] : enumerate(m_surfaces.internal)) {
    sensitiveId.setSensitive(i + 1);
    s->assignGeometryId(sensitiveId);
  }

  // Check if it is a container or detector volume
  if (not m_volumes.internal.empty()) {
    // Detection if any of the volume has sub surfaces
    bool detectorVolume = false;
    for (auto v : volumes()) {
      // Would in principle qualify
      if (not v->surfaces().empty()) {
        detectorVolume = true;
        break;
      }
    }
    // Cross-check if no container is present
    for (auto v : volumes()) {
      // Pure detector volume is vetoed
      if (not v->volumes().empty()) {
        detectorVolume = false;
        break;
      }
    }

    // Assign the volume Identifier (recursive step down)
    for (auto [i, v] : enumerate(m_volumes.internal)) {
      GeometryIdentifier volumeId = geometryId;
      if (detectorVolume) {
        volumeId.setLayer(i + 1);
      } else {
        volumeId.setVolume(volumeId.volume() + i + 1);
      }
      v->lock(volumeId);
    }
  }
}
