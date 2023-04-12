// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/DetectorVolume.hpp"

#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cassert>

Acts::Experimental::DetectorVolume::DetectorVolume(
    const GeometryContext& gctx, const std::string& name,
    const Transform3& transform, std::unique_ptr<VolumeBounds> bounds,
    const std::vector<std::shared_ptr<Surface>>& surfaces,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    SurfaceCandidatesUpdator&& surfaceCandidateUpdator)
    : m_name(name),
      m_transform(transform),
      m_bounds(std::move(bounds)),
      m_surfaceCandidatesUpdator(std::move(surfaceCandidateUpdator)),
      m_volumeMaterial(nullptr) {
  if (m_bounds == nullptr) {
    throw std::invalid_argument(
        "DetectorVolume: construction with nullptr bounds.");
  }
  if (not m_surfaceCandidatesUpdator.connected()) {
    throw std::invalid_argument(
        "DetectorVolume: navigation state updator delegate is not connected.");
  }

  m_surfaces = ObjectStore<std::shared_ptr<Surface>>(surfaces);
  m_volumes = ObjectStore<std::shared_ptr<DetectorVolume>>(volumes);

  [[maybe_unused]] const auto& gctx_ref = gctx;
  assert(checkContainment(gctx) && "Objects are not contained by volume.");
}

Acts::Experimental::DetectorVolume::DetectorVolume(
    const GeometryContext& /*gctx*/, const std::string& name,
    const Transform3& transform, std::unique_ptr<VolumeBounds> bounds,
    SurfaceCandidatesUpdator&& surfaceCandidateUpdator)
    : m_name(name),
      m_transform(transform),
      m_bounds(std::move(bounds)),
      m_surfaceCandidatesUpdator(std::move(surfaceCandidateUpdator)),
      m_volumeMaterial(nullptr) {
  if (m_bounds == nullptr) {
    throw std::invalid_argument(
        "DetectorVolume: construction with nullptr bounds.");
  }
  if (not m_surfaceCandidatesUpdator.connected()) {
    throw std::invalid_argument(
        "DetectorVolume: navigation state updator delegate is not connected.");
  }
}

void Acts::Experimental::DetectorVolume::updatePortal(
    std::shared_ptr<Portal> portal, unsigned int pIndex) {
  if (pIndex >= m_portals.internal.size()) {
    throw std::invalid_argument(
        "DetectorVolume: trying to update a portal that does not exist.");
  }
  m_portals.internal[pIndex] = std::move(portal);
  m_portals = ObjectStore<std::shared_ptr<Portal>>(m_portals.internal);
}

void Acts::Experimental::DetectorVolume::construct(
    const GeometryContext& gctx, const PortalGenerator& portalGenerator) {
  // Create portals with the given generator
  auto portalSurfaces =
      portalGenerator(transform(gctx), *(m_bounds.get()), getSharedPtr());
  m_portals = ObjectStore<std::shared_ptr<Portal>>(portalSurfaces);
  createBoundingBox(gctx);
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
                                                bool excludeInserts) const {
  Vector3 posInVolFrame((transform(gctx).inverse()) * position);
  if (not volumeBounds().inside(posInVolFrame)) {
    return false;
  }
  if (not excludeInserts or m_volumes.external.empty()) {
    return true;
  }
  // Check exclusion through subvolume
  for (const auto v : volumes()) {
    if (v->inside(gctx, position)) {
      return false;
    }
  }
  return true;
}

void Acts::Experimental::DetectorVolume::updateNavigationState(
    const GeometryContext& gctx, NavigationState& nState) const {
  nState.currentVolume = this;
  m_surfaceCandidatesUpdator(gctx, nState);
  nState.surfaceCandidate = nState.surfaceCandidates.begin();
}

void Acts::Experimental::DetectorVolume::assignSurfaceCandidatesUpdator(
    SurfaceCandidatesUpdator&& surfaceCandidateUpdator,
    const std::vector<std::shared_ptr<Surface>>& surfaces,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes) {
  m_surfaceCandidatesUpdator = std::move(surfaceCandidateUpdator);
  m_surfaces = ObjectStore<std::shared_ptr<Surface>>(surfaces);
  m_volumes = ObjectStore<std::shared_ptr<DetectorVolume>>(volumes);
}

Acts::Extent Acts::Experimental::DetectorVolume::extent(
    const GeometryContext& gctx, size_t nseg) const {
  Extent volumeExtent;
  for (const auto* p : portals()) {
    volumeExtent.extend(
        p->surface().polyhedronRepresentation(gctx, nseg).extent());
  }
  return volumeExtent;
}

bool Acts::Experimental::DetectorVolume::checkContainment(
    const GeometryContext& gctx, size_t nseg) const {
  // Create the volume extent
  auto volumeExtent = extent(gctx, nseg);
  // Check surfaces
  for (const auto* s : surfaces()) {
    auto sExtent = s->polyhedronRepresentation(gctx, nseg).extent();
    if (not volumeExtent.contains(sExtent)) {
      return false;
    }
  }
  // Check volumes
  for (const auto* v : volumes()) {
    auto vExtent = v->extent(gctx, nseg);
    if (not volumeExtent.contains(vExtent)) {
      return false;
    }
  }
  // All contained
  return true;
}

void Acts::Experimental::DetectorVolume::closePortals() {
  for (auto& p : m_portals.internal) {
    // Create a null link
    for (auto [ivu, vu] : enumerate(p->detectorVolumeUpdators())) {
      if (not vu.connected()) {
        auto eowDir = Acts::directionFromIndex(ivu);
        auto eow = std::make_unique<const EndOfWorldImpl>();
        Acts::Experimental::DetectorVolumeUpdator eowLink;
        eowLink.connect<&EndOfWorldImpl::update>(std::move(eow));
        p->assignDetectorVolumeUpdator(eowDir, std::move(eowLink), {});
      }
    }
  }

  for (auto& v : m_volumes.internal) {
    v->closePortals();
  }
}

void Acts::Experimental::DetectorVolume::createBoundingBox(
    const GeometryContext& gctx) {
  std::vector<Vector3> vertices;
  for (auto p : m_portals.external) {
    auto surface = p->surface().polyhedronRepresentation(gctx, 1);
    auto pVertices = surface.vertices;
    for (const auto& v : pVertices) {
      vertices.push_back(v);
    }
  }
  Acts::Vector3 vmin = Acts::Vector3::Zero();
  Acts::Vector3 vmax = Acts::Vector3::Zero();
  for (const auto& v : vertices) {
    vmin = vmin.cwiseMin(v);
    vmax = vmax.cwiseMax(v);
  }
  std::shared_ptr<Acts::Experimental::DetectorVolume::BoundingBox> box =
      std::make_shared<Acts::Experimental::DetectorVolume::BoundingBox>(
          this, vmin, vmax);
  m_boundingBox = box;
}

const Acts::Experimental::DetectorVolume::BoundingBox&
Acts::Experimental::DetectorVolume::getBoundingBox() const {
  return *m_boundingBox;
}
