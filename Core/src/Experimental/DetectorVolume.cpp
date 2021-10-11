// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/DetectorVolume.hpp"

#include "Acts/Experimental/Enumerate.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/PortalLinks.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"

Acts::DetectorVolume::DetectorVolume(const Transform3& transform,
                                     std::unique_ptr<VolumeBounds> bounds,
                                     const std::string& name)
    : m_transform(transform), m_bounds(std::move(bounds)), m_name(name) {
  // This will create the portal surfaces
  createPortals();
}

Acts::DetectorVolume::DetectorVolume(
    const Transform3& transform, std::unique_ptr<VolumeBounds> bounds,
    std::vector<std::shared_ptr<Surface>>&& surfaces,
    SurfaceLinks&& volumeSurfaceLinks,
    std::vector<SurfaceLinks>&& portalSurfaceLinks, const std::string& name)
    : m_transform(transform),
      m_bounds(std::move(bounds)),
      m_surfaces(std::move(surfaces)),
      m_surfaceLinks(volumeSurfaceLinks),
      m_name(name) {
  // This will create the portal surfaces
  createPortals(std::move(portalSurfaceLinks));
}

Acts::DetectorVolume::DetectorVolume(
    const Transform3& transform, std::unique_ptr<VolumeBounds> bounds,
    std::vector<std::shared_ptr<DetectorVolume>>&& volumes,
    VolumeLink&& volumeLink, bool recastPortals, BinningValue recastValue,
    const std::string& name)
    : m_transform(transform),
      m_bounds(std::move(bounds)),
      m_volumes(std::move(volumes)),
      m_volumeLink(volumeLink),
      m_name(name) {
  // This will create the portal surfaces
  createPortals({}, recastPortals, recastValue);
}

bool Acts::DetectorVolume::inside(const Vector3& position,
                                  ActsScalar tolerance) const {
  Vector3 posInVolFrame((transform().inverse()) * position);
  return (volumeBounds()).inside(posInVolFrame, tolerance);
}

void Acts::DetectorVolume::createPortals(
    std::vector<SurfaceLinks>&& surfaceLinks, bool recastPortals,
    BinningValue recastValue) {
  auto orientedSurfaces = m_bounds->orientedSurfaces(m_transform);
  // Create & reserve portal vector
  std::vector<std::shared_ptr<Portal>> portalSurfaces;
  portalSurfaces.reserve(orientedSurfaces.size());
  // No recasting of portals for a container to be done
  if (not recastPortals) {
    // Create surfaces and translate them into portals
    if (not surfaceLinks.empty() and
        orientedSurfaces.size() != surfaceLinks.size()) {
      throw std::invalid_argument(
          "DetectorVolume: not enough surface link objects provided");
    }

    // Loop over the oriented surfaces vector and
    for (auto [i, osf] : enumerate(orientedSurfaces)) {
      // Set the volume to the surface link
      SinglePortalLink pLink = SinglePortalLink{this};
      if (not surfaceLinks.empty()) {
        pLink.surfaces = std::move(surfaceLinks[i]);
      }
      auto portalSurface = std::make_shared<Portal>(osf.first);
      portalSurface->updatePortalLink(std::move(pLink), osf.second);
      portalSurfaces.push_back(std::move(portalSurface));
    }
  } else {
    // Recasting option chosen
    if (m_bounds->type() == VolumeBounds::eCylinder) {
      size_t nsf = orientedSurfaces.size();
      // The internal volume representation
      auto ivolumes = m_volumes.internal;
      auto volume0 = ivolumes[0];
      auto volumeL = ivolumes[ivolumes.size() - 1];
      if (recastValue == binR) {
        portalSurfaces.push_back(volume0->m_portals.internal[0]);
        portalSurfaces.push_back(volume0->m_portals.internal[1]);
        portalSurfaces.push_back(volumeL->m_portals.internal[2]);
        if (nsf > 3u) {
          portalSurfaces.push_back(volume0->m_portals.internal[3]);
        }
        // @TODO include phi sectors
      } else if (recastValue == binZ) {
        portalSurfaces.push_back(volume0->m_portals.internal[0]);
        portalSurfaces.push_back(volumeL->m_portals.internal[1]);
        portalSurfaces.push_back(volume0->m_portals.internal[2]);
        if (nsf > 3u) {
          portalSurfaces.push_back(volume0->m_portals.internal[3]);
        }
        // @TODO include phi sectors
      }
    } else {
      throw std::invalid_argument("DetectorVolume: not yet implemented.");
    }
  }
  m_portals = ObjectStore<std::shared_ptr<Portal>>(std::move(portalSurfaces));
}

void Acts::DetectorVolume::updatePortalPtr(
    std::shared_ptr<Portal> updatedPortal, size_t portal, bool keepPortalLink,
    NavigationDirection nDir) {
  // Throw an exception if no portals are present
  if (portal >= m_portals.internal.size()) {
    throw std::invalid_argument(
        "DetectorVolume: trying to update non-existing portal.");
  }

  // If configured, keep the portal link, DetectorVolume is
  // friend of the Portal, so it is allowed to do so
  if (keepPortalLink and nDir == forward) {
    updatedPortal->m_alongNormal =
        std::move(m_portals.internal[portal]->m_alongNormal);
  } else if (keepPortalLink) {
    updatedPortal->m_oppositeNormal =
        std::move(m_portals.internal[portal]->m_oppositeNormal);
  }
  // Now overwrite the the portal
  m_portals.internal[portal] = updatedPortal;
  // Reninitalize the portal store
  m_portals =
      ObjectStore<std::shared_ptr<Portal>>(std::move(m_portals.internal));
}

Acts::DetectorEnvironment Acts::DetectorVolume::environment(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bCheck,
    bool provideAll) const {
  DetectorEnvironment nEnvironment{};

  // The volume is obviously this
  nEnvironment.volume = this;
  // Evaluate the portals
  nEnvironment.portals = portalCandidates(gctx, *this, position, direction);

  // Maximum path length, from closest portal intersection
  ActsScalar maximumPath = nEnvironment.portals[0].intersection.pathLength;

  // The surface candidates one-time intersected & ordered
  nEnvironment.surfaces = m_surfaceLinks(gctx, *this, position, direction,
                                         bCheck, maximumPath, provideAll);

  return nEnvironment;
}

const Acts::DetectorVolume* Acts::DetectorVolume::lowest(
    const GeometryContext& gctx, const Vector3& position) const {
  const auto& vs = volumes();
  if (not vs.empty()) {
    unsigned int v = m_volumeLink(transform(gctx), position);
    return vs[v]->lowest(gctx, position);
  }
  return this;
}
