// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/DetectorVolume.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <cassert>
#include <iterator>

namespace Acts {
class IVolumeMaterial;
}  // namespace Acts

Acts::Experimental::DetectorVolume::DetectorVolume(
    const GeometryContext& gctx, std::string name, const Transform3& transform,
    std::shared_ptr<VolumeBounds> bounds,
    std::vector<std::shared_ptr<Surface>> surfaces,
    std::vector<std::shared_ptr<DetectorVolume>> volumes,
    DetectorVolumeUpdator detectorVolumeUpdator,
    SurfaceCandidatesUpdator surfaceCandidateUpdator)
    : m_name(std::move(name)),
      m_transform(transform),
      m_bounds(std::move(bounds)),
      m_surfaces(std::move(surfaces)),
      m_volumes(std::move(volumes)),
      m_detectorVolumeUpdator(std::move(detectorVolumeUpdator)),
      m_surfaceCandidatesUpdator(std::move(surfaceCandidateUpdator)),
      m_volumeMaterial(nullptr) {
  if (m_bounds == nullptr) {
    throw std::invalid_argument(
        "DetectorVolume: construction with nullptr bounds.");
  }
  if (not m_detectorVolumeUpdator.connected()) {
    throw std::invalid_argument(
        "DetectorVolume: navigation state updator delegate is not connected.");
  }
  if (not m_surfaceCandidatesUpdator.connected()) {
    throw std::invalid_argument(
        "DetectorVolume: navigation state updator delegate is not connected.");
  }

  [[maybe_unused]] const auto& gctx_ref = gctx;
  assert(checkContainment(gctx) && "Objects are not contained by volume.");
}

Acts::Experimental::DetectorVolume::DetectorVolume(
    const GeometryContext& gctx, std::string name, const Transform3& transform,
    std::shared_ptr<VolumeBounds> bounds,
    SurfaceCandidatesUpdator surfaceCandidateUpdator)
    : DetectorVolume(gctx, std::move(name), transform, std::move(bounds), {},
                     {}, tryNoVolumes(), std::move(surfaceCandidateUpdator)) {}

std::shared_ptr<Acts::Experimental::DetectorVolume>
Acts::Experimental::DetectorVolume::makeShared(
    const GeometryContext& gctx, std::string name, const Transform3& transform,
    std::shared_ptr<VolumeBounds> bounds,
    std::vector<std::shared_ptr<Surface>> surfaces,
    std::vector<std::shared_ptr<DetectorVolume>> volumes,
    DetectorVolumeUpdator detectorVolumeUpdator,
    SurfaceCandidatesUpdator surfaceCandidateUpdator) {
  return std::shared_ptr<DetectorVolume>(new DetectorVolume(
      gctx, std::move(name), transform, std::move(bounds), std::move(surfaces),
      std::move(volumes), std::move(detectorVolumeUpdator),
      std::move(surfaceCandidateUpdator)));
}

std::shared_ptr<Acts::Experimental::DetectorVolume>
Acts::Experimental::DetectorVolume::makeShared(
    const GeometryContext& gctx, std::string name, const Transform3& transform,
    std::shared_ptr<VolumeBounds> bounds,
    SurfaceCandidatesUpdator surfaceCandidateUpdator) {
  return std::shared_ptr<DetectorVolume>(
      new DetectorVolume(gctx, std::move(name), transform, std::move(bounds),
                         std::move(surfaceCandidateUpdator)));
}

const Acts::Transform3& Acts::Experimental::DetectorVolume::transform(
    const GeometryContext& /*gctx*/) const {
  return m_transform;
}

Acts::Vector3 Acts::Experimental::DetectorVolume::center(
    const GeometryContext& gctx) const {
  return transform(gctx).translation();
}

const Acts::VolumeBounds& Acts::Experimental::DetectorVolume::volumeBounds()
    const {
  return (*m_bounds.get());
}

std::vector<std::shared_ptr<Acts::Experimental::Portal>>&
Acts::Experimental::DetectorVolume::portalPtrs() {
  return m_portals.internal;
}

std::vector<std::shared_ptr<Acts::Surface>>&
Acts::Experimental::DetectorVolume::surfacePtrs() {
  return m_surfaces.internal;
}

std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>&
Acts::Experimental::DetectorVolume::volumePtrs() {
  return m_volumes.internal;
}

const std::vector<const Acts::Experimental::Portal*>&
Acts::Experimental::DetectorVolume::portals() const {
  return m_portals.external;
}

const std::vector<const Acts::Surface*>&
Acts::Experimental::DetectorVolume::surfaces() const {
  return m_surfaces.external;
}

const std::vector<const Acts::Experimental::DetectorVolume*>&
Acts::Experimental::DetectorVolume::volumes() const {
  return m_volumes.external;
}

const Acts::Experimental::DetectorVolumeUpdator&
Acts::Experimental::DetectorVolume::detectorVolumeUpdator() const {
  return m_detectorVolumeUpdator;
}

const Acts::Experimental::SurfaceCandidatesUpdator&
Acts::Experimental::DetectorVolume::surfaceCandidatesUpdator() const {
  return m_surfaceCandidatesUpdator;
}

void Acts::Experimental::DetectorVolume::assignVolumeMaterial(
    std::shared_ptr<IVolumeMaterial> material) {
  m_volumeMaterial = std::move(material);
}

std::shared_ptr<Acts::IVolumeMaterial>
Acts::Experimental::DetectorVolume::volumeMaterialPtr() {
  return m_volumeMaterial;
}

const Acts::IVolumeMaterial*
Acts::Experimental::DetectorVolume::volumeMaterial() const {
  return m_volumeMaterial.get();
}

const Acts::GeometryIdentifier& Acts::Experimental::DetectorVolume::geometryId()
    const {
  return m_geometryId;
}

void Acts::Experimental::DetectorVolume::assignGeometryId(
    const GeometryIdentifier& geoID) {
  m_geometryId = geoID;
}

const std::string& Acts::Experimental::DetectorVolume::name() const {
  return m_name;
}

void Acts::Experimental::DetectorVolume::assignDetector(
    const Detector& detector) {
  m_detector = &detector;

  for (auto& v : m_volumes.internal) {
    v->assignDetector(detector);
  }
}

const Acts::Experimental::Detector*
Acts::Experimental::DetectorVolume::detector() const {
  return m_detector;
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
                                                const Vector3& position) const {
  Vector3 posInVolFrame(transform(gctx).inverse() * position);
  return volumeBounds().inside(posInVolFrame);
}

bool Acts::Experimental::DetectorVolume::exclusivelyInside(
    const GeometryContext& gctx, const Vector3& position) const {
  if (!inside(gctx, position)) {
    return false;
  }
  // Check exclusion through subvolume
  for (const auto& v : volumes()) {
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
    SurfaceCandidatesUpdator surfaceCandidateUpdator,
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
        auto eowDir = Direction::fromIndex(ivu);
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
