// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylindricalDetectorHelper.hpp"

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Detector.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Geometry/detail/DetectorVolumeFinders.hpp"
#include "Acts/Geometry/detail/PortalHelper.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <exception>
#include <unordered_set>

namespace {

/// @brief Helper function to create a disc portal replacement
/// @param transform The transform of the newly created portal
/// @param rBoundaris the vector of binning boundaries in r
/// @param phiSector an eventual phi sector value
/// @param avgPhi an eventual average phi parameter
/// @param index the index of this portal to be set
/// @param nDir the navigation direction to be set
///
/// @return a new portal replacement object
Acts::Experimental::detail::PortalReplacement createDiscReplacement(
    const Acts::Transform3& transform,
    const std::vector<Acts::ActsScalar>& rBoundaries,
    Acts::ActsScalar phiSector, Acts::ActsScalar avgPhi, unsigned int index,
    Acts::NavigationDirection nDir) {
  // Transform and bounds
  auto bounds = std::make_unique<Acts::RadialBounds>(
      rBoundaries[0u], rBoundaries[rBoundaries.size() - 1u], phiSector, avgPhi);
  // A new surface on the negative side over the full range
  auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(
      transform, std::move(bounds));
  // Make a portal and indicate the new link direction
  Acts::Experimental::detail::PortalReplacement pRep = {
      Acts::Experimental::Portal::makeShared(surface), index, nDir, rBoundaries,
      Acts::binR};
  return pRep;
}

/// @brief Helper function to create a disc portal replacement
/// @param transform The transform of the newly created portal
/// @param boundaries the vector of binning boundaries
/// @param r is the radius of the cylinder
/// @param phiSector an eventual phi sector value
/// @param avgPhi an eventual average phi parameter
/// @param index the index of this portal to be set
/// @param nDir the navigation direction to be set
///
/// @return a new portal replacement object
Acts::Experimental::detail::PortalReplacement createCylinderReplacement(
    const Acts::Transform3& transform,
    const std::vector<Acts::ActsScalar>& boundaries, Acts::ActsScalar r,
    Acts::ActsScalar phiSector, Acts::ActsScalar avgPhi, unsigned int index,
    Acts::NavigationDirection nDir) {
  Acts::ActsScalar diffZ = boundaries[boundaries.size() - 1u] - boundaries[0u];
  // Transform and bounds
  auto bounds =
      std::make_unique<Acts::CylinderBounds>(r, 0.5 * diffZ, phiSector, avgPhi);
  // Re-scale the boundaries
  ///@todo needs to change to an appropriate transform
  // A new surface on the negative side over the full range
  auto surface = Acts::Surface::makeShared<Acts::CylinderSurface>(
      transform, std::move(bounds));
  // A make a portal and indicate the new link direction
  Acts::Experimental::detail::PortalReplacement pRep = {
      Acts::Experimental::Portal::makeShared(surface), index, nDir, boundaries,
      Acts::binZ};

  return pRep;
}

/// @brief Helper function to create a disc portal replacement
/// @param gcxt the geometry context of thix call
/// @param volumeCenter a reference center of the volume (combined)
/// @param refSurface the reference surface (old portal)
/// @param boundaries the vector of binning boundaries in r
/// @param binning the binning of the sector (inR, inZ)
/// @param index the index of this portal to be set
/// @param nDir the navigation direction to be set
///
/// @return a new portal replacement object
Acts::Experimental::detail::PortalReplacement createSectorReplacement(
    const Acts::GeometryContext& gctx, const Acts::Vector3& volumeCenter,
    const Acts::Surface& refSurface,
    const std::vector<Acts::ActsScalar>& boundaries, Acts::BinningValue binning,
    unsigned int index, Acts::NavigationDirection nDir) {
  // Get a reference transform
  const auto& refTransform = refSurface.transform(gctx);
  auto refRotation = refTransform.rotation();
  // Bounds handling
  const auto& boundValues = refSurface.bounds().values();
  std::unique_ptr<Acts::PlanarBounds> bounds = nullptr;

  // Create a new transform
  Acts::Transform3 transform = Acts::Transform3::Identity();
  if (binning == Acts::binR) {
    // Center r calculation
    Acts::ActsScalar centerR =
        0.5 * (boundaries[0u] + boundaries[boundaries.size() - 1u]);
    // New joint center
    Acts::Vector3 pCenter = volumeCenter + centerR * refRotation.col(2);
    transform.pretranslate(pCenter);
    Acts::ActsScalar halfX =
        0.5 * (boundValues[Acts::RectangleBounds::BoundValues::eMaxX] -
               boundValues[Acts::RectangleBounds::BoundValues::eMinX]);
    // New joint bounds
    bounds = std::make_unique<Acts::RectangleBounds>(
        halfX, 0.5 * (boundaries[boundaries.size() - 1u] - boundaries[0u]));
  } else if (binning == Acts::binZ) {
    // Center R calculation, using projection onto vector
    const auto& surfaceCenter = refSurface.center(gctx);
    Acts::Vector3 centerDiffs = (surfaceCenter - volumeCenter);
    Acts::ActsScalar centerR = centerDiffs.dot(refRotation.col(2));
    // New joint center
    Acts::Vector3 pCenter = volumeCenter + centerR * refRotation.col(2);
    transform.pretranslate(pCenter);
    // New joint bounds
    Acts::ActsScalar halfY =
        0.5 * (boundValues[Acts::RectangleBounds::BoundValues::eMaxY] -
               boundValues[Acts::RectangleBounds::BoundValues::eMinY]);
    bounds = std::make_unique<Acts::RectangleBounds>(
        0.5 * (boundaries[boundaries.size() - 1u] - boundaries[0u]), halfY);
  }
  // Set the rotation
  transform.prerotate(refRotation);
  // A new surface on the negative side over the full range
  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform, std::move(bounds));
  // A make a portal and indicate the new link direction
  Acts::Experimental::detail::PortalReplacement pRep = {
      Acts::Experimental::Portal::makeShared(surface), index, nDir, boundaries,
      binning};
  return pRep;
}

/// @brief  Helper method to strip side volumes from containers
///
/// @param containers the list of container
/// @param sides the sides
/// @param selectedOnly the selection restriction
///
/// @return a map of stripped out container
std::map<unsigned int,
         std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>>
stripSideVolumes(
    const std::vector<Acts::Experimental::ProtoContainer>& containers,
    const std::vector<unsigned int>& sides,
    const std::vector<unsigned int>& selectedOnly = {}) {
  // Thes are the aripped out side volumes
  std::map<unsigned int,
           std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>>
      sideVolumes;

  // Principle sides and selected sides, make an intersection
  std::vector<unsigned int> selectedSides;
  if (not selectedOnly.empty()) {
    std::set_intersection(sides.begin(), sides.end(), selectedOnly.begin(),
                          selectedOnly.end(),
                          std::back_inserter(selectedSides));
  } else {
    selectedSides = sides;
  }

  for (const auto& pc : containers) {
    for (const auto& s : selectedSides) {
      if (pc.find(s) != pc.end()) {
        auto p = pc.find(s)->second;
        if (sideVolumes.find(s) == sideVolumes.end()) {
          sideVolumes[s] = {};
        }
        auto& sVolumes = sideVolumes[s];
        auto aVolumes =
            Acts::Experimental::detail::attachedDetectorVolumes(*p.get());
        sVolumes.insert(sVolumes.end(), aVolumes.begin(), aVolumes.end());
      }
    }
  }
  // return them
  return sideVolumes;
}
}  // namespace

Acts::Experimental::ProtoContainer
Acts::Experimental::connectDetectorVolumesInR(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<Experimental::DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << volumes.size() << " detector volumes in R.");

  // Return proto container
  ProtoContainer protoContainer;

  // Innermost volume boundaries
  std::vector<ActsScalar> rBoundaries = {};
  auto refValues = volumes[0u]->volumeBounds().values();

  // Reference boundary values
  rBoundaries.push_back(refValues[0u]);
  rBoundaries.push_back(refValues[1u]);

  // Connect in R ?
  bool connectR = selectedOnly.empty() or
                  std::find(selectedOnly.begin(), selectedOnly.end(), 2u) !=
                      selectedOnly.end();

  // Get phi sector and average phi
  ActsScalar phiSector =
      refValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector];
  ActsScalar avgPhi = refValues[CylinderVolumeBounds::BoundValues::eAveragePhi];

  // Fuse the cylinders - portals can be reused for this operation
  for (unsigned int iv = 1; iv < volumes.size(); ++iv) {
    // collect the radial bounds
    refValues = volumes[iv]->volumeBounds().values();
    rBoundaries.push_back(refValues[1u]);
    // Only connect if configured to do so
    if (connectR) {
      // Fuse and swap
      auto& keepCylinder = volumes[iv - 1]->portalPtrs()[2u];
      auto& wasteCylinder = volumes[iv]->portalPtrs()[3u];
      keepCylinder->fuse(wasteCylinder);
      volumes[iv]->updatePortal(keepCylinder, 3u);
    }
  }

  // Proto container setting
  if (connectR) {
    if (volumes[0u]->portals().size() == 4u or
        volumes[0u]->portals().size() == 6u) {
      protoContainer[3u] = volumes[0u]->portalPtrs()[3u];
    }
    protoContainer[2u] = volumes[volumes.size() - 1u]->portalPtrs()[2u];
  }

  // Check if sectors are present by the number of portals
  bool sectorsPresent = volumes[volumes.size() - 1u]->portals().size() > 4u;

  // A portal replacement, it comprises of the portal, the index, the direction,
  // the binning and bins
  std::vector<detail::PortalReplacement> pReplacements = {};

  // Disc assignments are forwad for negative disc, backward for positive
  std::vector<Acts::NavigationDirection> discDirs = {
      Acts::NavigationDirection::Forward, Acts::NavigationDirection::Backward};
  for (const auto [iu, idir] : enumerate(discDirs)) {
    if (selectedOnly.empty() or
        std::find(selectedOnly.begin(), selectedOnly.end(), iu) !=
            selectedOnly.end()) {
      const Surface& refSurface = volumes[0u]->portals()[iu]->surface();
      const Transform3& refTransform = refSurface.transform(gctx);
      pReplacements.push_back(createDiscReplacement(
          refTransform, rBoundaries, phiSector, avgPhi, iu, idir));
    }
  }

  // If volume sectors are present, these have to be respected
  if (sectorsPresent) {
    // Sector assignmenta are forward backward
    std::vector<Acts::NavigationDirection> sectorDirs = {
        Acts::NavigationDirection::Forward,
        Acts::NavigationDirection::Backward};
    Acts::Vector3 vCenter = volumes[0u]->transform(gctx).translation();
    for (const auto [iu, idir] : enumerate(sectorDirs)) {
      if (selectedOnly.empty() or
          std::find(selectedOnly.begin(), selectedOnly.end(), iu + 4u) !=
              selectedOnly.end()) {
        // As it is r-wrapping, the inner tube is guaranteed
        const Surface& refSurface =
            volumes[volumes.size() - 1u]->portals()[iu + 4u]->surface();
        pReplacements.push_back(createSectorReplacement(
            gctx, vCenter, refSurface, rBoundaries, Acts::binR, iu + 4u, idir));
      }
    }
  }

  // Attach the new volume multi links
  detail::attachDetectorVolumeUpdators(gctx, volumes, pReplacements);

  // Exchange the portals of the volumes
  for (auto& iv : volumes) {
    for (auto& [p, i, dir, boundaries, binning] : pReplacements) {
      // Fill the map
      protoContainer[i] = p;
      // Potential offset for tube vs/ cylinder
      int iOffset = (iv->portals().size() == 3u and i > 2u) ? -1 : 0;
      iv->updatePortal(p, static_cast<unsigned int>(i + iOffset));
    }
  }
  // Done.
  return protoContainer;
}

Acts::Experimental::ProtoContainer
Acts::Experimental::connectDetectorVolumesInZ(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << volumes.size() << " detector volumes in z.");

  // Return proto container
  ProtoContainer protoContainer;

  // Connect in Z ?
  bool connectZ = selectedOnly.empty() or
                  std::find(selectedOnly.begin(), selectedOnly.end(), 1u) !=
                      selectedOnly.end();
  // Reference z axis
  const auto rotation = volumes[0]->transform(gctx).rotation();
  std::vector<Vector3> zBoundaries3D = {};

  /// @brief  Add the z boundary
  /// @param gctx the geometry context
  /// @param volume the volume
  /// @param side side value
  auto addZboundary3D = [&](const Experimental::DetectorVolume& volume,
                            int side) -> void {
    const auto boundValues = volume.volumeBounds().values();
    ActsScalar halflengthZ =
        boundValues[CylinderVolumeBounds::BoundValues::eHalfLengthZ];
    zBoundaries3D.push_back(volume.transform(gctx).translation() +
                            side * halflengthZ * rotation.col(2));
  };

  // Fuse the discs - portals can be reused
  addZboundary3D(*volumes[0u].get(), -1);
  addZboundary3D(*volumes[0u].get(), 1);
  for (unsigned int iv = 1; iv < volumes.size(); ++iv) {
    // Add the z boundary
    addZboundary3D(*volumes[iv].get(), 1u);
    // Do the connection
    if (connectZ) {
      auto& keepDisc = volumes[iv - 1]->portalPtrs()[1u];
      auto& wasteDisc = volumes[iv]->portalPtrs()[0u];
      keepDisc->fuse(wasteDisc);
      volumes[iv]->updatePortal(keepDisc, 0u);
    }
  }

  // Register what we have from the container
  if (connectZ) {
    protoContainer[0u] = volumes[0u]->portalPtrs()[0u];
    protoContainer[1u] = volumes[volumes.size() - 1u]->portalPtrs()[1u];
  }

  // Centre of gravity
  Vector3 combinedCenter =
      0.5 * (zBoundaries3D[zBoundaries3D.size() - 1u] + zBoundaries3D[0u]);

  std::vector<ActsScalar> zBoundaries = {};
  for (const auto& zb3D : zBoundaries3D) {
    auto proj3D = (zb3D - combinedCenter).dot(rotation.col(2));
    ActsScalar zBoundary =
        std::copysign((zb3D - combinedCenter).norm(), proj3D);
    zBoundaries.push_back(zBoundary);
  }

  Transform3 combinedTransform = Transform3::Identity();
  combinedTransform.pretranslate(combinedCenter);
  combinedTransform.rotate(rotation);

  // Get phi sector and average phi
  const auto& refVolume = volumes[0u];
  const auto refValues = refVolume->volumeBounds().values();
  // Get phi sector and average phi
  ActsScalar minR = refValues[CylinderVolumeBounds::BoundValues::eMinR];
  ActsScalar maxR = refValues[CylinderVolumeBounds::BoundValues::eMaxR];
  ActsScalar phiSector =
      refValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector];
  ActsScalar avgPhi = refValues[CylinderVolumeBounds::BoundValues::eAveragePhi];

  // Check if inner cylinder and sectors are present by the number of portals
  size_t nPortals = volumes[volumes.size() - 1u]->portals().size();
  bool innerPresent = (nPortals != 3u and nPortals != 5u);
  bool sectorsPresent = nPortals > 4u;

  // A portal replacement, it comprises of the portal, the index, the direction,
  // the binning and bins
  std::vector<detail::PortalReplacement> pReplacements = {};

  // Disc assignments are forwad for negative disc, backward for positive
  std::vector<Acts::NavigationDirection> cylinderDirs = {
      Acts::NavigationDirection::Backward};
  // Cylinder radii
  std::vector<Acts::ActsScalar> cylinderR = {maxR};
  if (innerPresent) {
    cylinderDirs.push_back(NavigationDirection::Forward);
    cylinderR.push_back(minR);
  }
  // Tube/cylinder offset
  unsigned int iSecOffset = innerPresent ? 4u : 3u;
  // Prepare the cylinder replacements
  for (const auto [iu, idir] : enumerate(cylinderDirs)) {
    if (selectedOnly.empty() or
        std::find(selectedOnly.begin(), selectedOnly.end(), iu + 2u) !=
            selectedOnly.end()) {
      pReplacements.push_back(createCylinderReplacement(
          combinedTransform, zBoundaries, cylinderR[iu], phiSector, avgPhi,
          iu + 2u, idir));
    }
  }

  // Prepare the sector side replacements
  if (sectorsPresent) {
    // Sector assignmenta are forward backward
    std::vector<Acts::NavigationDirection> sectorDirs = {
        Acts::NavigationDirection::Forward,
        Acts::NavigationDirection::Backward};
    for (const auto [iu, idir] : enumerate(sectorDirs)) {
      // Access with 3u or 4u but always write 4u (to be caught later)
      if (selectedOnly.empty() or
          std::find(selectedOnly.begin(), selectedOnly.end(), iu + 4u) !=
              selectedOnly.end()) {
        const Surface& refSurface =
            volumes[0u]->portals()[iu + iSecOffset]->surface();
        pReplacements.push_back(
            createSectorReplacement(gctx, combinedCenter, refSurface,
                                    zBoundaries, Acts::binZ, iu + 4u, idir));
      }
    }
  }

  // Attach the new volume multi links
  detail::attachDetectorVolumeUpdators(gctx, volumes, pReplacements);

  // Exchange the portals of the volumes
  for (auto& iv : volumes) {
    for (auto& [p, i, dir, boundaries, binning] : pReplacements) {
      // Potential offset for tube vs/ cylinder
      int iOffset = (i > 2u and not innerPresent) ? -1 : 0;
      iv->updatePortal(p, static_cast<unsigned int>(i + iOffset));
      // Fill the map
      protoContainer[i] = p;
    }
  }
  // Done.
  return protoContainer;
}

Acts::Experimental::ProtoContainer
Acts::Experimental::connectDetectorVolumesInPhi(
    const Acts::GeometryContext& /*gctx*/,
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>& volumes,
    const std::vector<unsigned int>& /*selectedOnly*/,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << volumes.size() << " detector volumes in phi.");

  // Return proto container
  ProtoContainer protoContainer;

  // Check if inner cylinder and sectors are present by the number of portals
  size_t nPortals = volumes[volumes.size() - 1u]->portals().size();
  bool innerPresent = (nPortals != 3u and nPortals != 5u);

  unsigned int iSecOffset = innerPresent ? 4u : 3u;
  // Fuse the cylinders - portals can be reused for this operation
  for (unsigned int iv = 1; iv < volumes.size(); ++iv) {
    // Fuse and swap
    auto& keepSector = volumes[iv - 1]->portalPtrs()[iSecOffset + 1u];
    auto& wasteSector = volumes[iv]->portalPtrs()[iSecOffset];
    keepSector->fuse(wasteSector);
    volumes[iv]->updatePortal(keepSector, iSecOffset);
  }
  // Done.
  return protoContainer;
}

Acts::Experimental::ProtoContainer Acts::Experimental::connectContainersInR(
    const GeometryContext& gctx, const std::vector<ProtoContainer>& containers,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) noexcept(false) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << containers.size() << " proto containers in R.");

  // Return the new container
  ProtoContainer protoContainer;

  // Fuse the cylinders - portals can be reused for this operation
  for (unsigned int ic = 1; ic < containers.size(); ++ic) {
    auto& formerContainer = containers[ic - 1];
    auto& currentContainer = containers[ic];
    // Check and throw exception
    if (formerContainer.find(2u) == formerContainer.end()) {
      throw std::invalid_argument(
          "CylindricalDetectorHelper: proto container has no outer cover, can "
          "not be connected in R");
    }
    if (currentContainer.find(3u) == currentContainer.end()) {
      throw std::invalid_argument(
          "CylindricalDetectorHelper: proto container has no inner cover, can "
          "not be connected in R");
    }

    // Fuse and swap
    std::shared_ptr<Portal> keepCylinder = containers[ic - 1].find(2u)->second;
    std::shared_ptr<Portal> wasteCylinder = containers[ic].find(3u)->second;
    keepCylinder->fuse(wasteCylinder);
    for (auto& av : wasteCylinder->attachedDetectorVolumes()[1u]) {
      av->updatePortal(keepCylinder, 3u);
    }
  }

  // Proto container refurbishment
  if (containers[0u].find(3u) != containers[0u].end()) {
    protoContainer[3u] = containers[0u].find(3u)->second;
  }
  protoContainer[2u] = containers[containers.size() - 1u].find(2u)->second;

  auto sideVolumes =
      stripSideVolumes(containers, {0u, 1u, 4u, 5u}, selectedOnly);

  for (auto [s, volumes] : sideVolumes) {
    auto pR = connectDetectorVolumesInR(gctx, volumes, {s});
    if (pR.find(s) != pR.end()) {
      protoContainer[s] = pR.find(s)->second;
    }
  }

  // Done.
  return protoContainer;
}

Acts::Experimental::ProtoContainer Acts::Experimental::connectContainersInZ(
    const GeometryContext& gctx, const std::vector<ProtoContainer>& containers,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) noexcept(false) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << containers.size() << " proto containers in Z.");

  // Return the new container
  ProtoContainer protoContainer;

  for (unsigned int ic = 1; ic < containers.size(); ++ic) {
    auto& formerContainer = containers[ic - 1];
    auto& currentContainer = containers[ic];
    // Check and throw exception
    if (formerContainer.find(1u) == formerContainer.end()) {
      throw std::invalid_argument(
          "CylindricalDetectorHelper: proto container has no negative disc, "
          "can "
          "not be connected in Z");
    }
    if (currentContainer.find(0u) == currentContainer.end()) {
      throw std::invalid_argument(
          "CylindricalDetectorHelper: proto container has no positive disc, "
          "can "
          "not be connected in Z");
    }
    std::shared_ptr<Portal> keepDisc = formerContainer.find(1u)->second;
    std::shared_ptr<Portal> wasteDisc = currentContainer.find(0u)->second;
    keepDisc->fuse(wasteDisc);
    for (auto& av : wasteDisc->attachedDetectorVolumes()[1u]) {
      av->updatePortal(keepDisc, 0u);
    }
  }

  // Proto container refurbishment
  protoContainer[0u] = containers[0u].find(0u)->second;
  protoContainer[1u] = containers[containers.size() - 1u].find(1u)->second;

  std::vector<unsigned int> nominalSides = {2u, 4u, 5u};
  if (containers[0u].find(3u) != containers[0u].end()) {
    nominalSides.push_back(3u);
  }

  // Set the sides
  auto sideVolumes = stripSideVolumes(containers, nominalSides, selectedOnly);

  for (auto [s, volumes] : sideVolumes) {
    auto pR = connectDetectorVolumesInZ(gctx, volumes, {s});
    if (pR.find(s) != pR.end()) {
      protoContainer[s] = pR.find(s)->second;
    }
  }

  // Done.
  return protoContainer;
}


std::array<std::vector<Acts::ActsScalar>, 3u>
Acts::Experimental::rzphiBoundaries(
    const GeometryContext& gctx,
    const std::vector<const Acts::Experimental::DetectorVolume*>& volumes,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Estimate R/z/phi boundaries of  " << volumes.size()
                                                << " volumes.");

  // The return boundaries
  std::array<std::vector<Acts::ActsScalar>, 3u> boundaries;

  // The map for collecting
  std::array<std::map<ActsScalar, size_t>, 3u> valueMaps;
  auto& rMap = valueMaps[0u];
  auto& zMap = valueMaps[1u];
  auto& phiMap = valueMaps[2u];

  auto fillMap = [&](std::map<ActsScalar, size_t>& map,
                     const std::array<ActsScalar, 2u>& values) {
    for (auto v : values) {
      if (map.find(v) != map.end()) {
        ++map[v];
      } else {
        map[v] = 1u;
      }
    }
  };

  // Loop over the volumes and collect boundaries
  for (const auto& v : volumes) {
    if (v->volumeBounds().type() == Acts::VolumeBounds::BoundsType::eCylinder) {
      auto bValues = v->volumeBounds().values();
      // The min/max values
      ActsScalar rMin = bValues[CylinderVolumeBounds::BoundValues::eMinR];
      ActsScalar rMax = bValues[CylinderVolumeBounds::BoundValues::eMaxR];
      ActsScalar zCenter = v->transform(gctx).translation().z();
      ActsScalar zHalfLength =
          bValues[CylinderVolumeBounds::BoundValues::eHalfLengthZ];
      ActsScalar zMin = zCenter - zHalfLength;
      ActsScalar zMax = zCenter + zHalfLength;
      ActsScalar phiCenter =
          bValues[CylinderVolumeBounds::BoundValues::eAveragePhi];
      ActsScalar phiSector =
          bValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector];
      ActsScalar phiMin = phiCenter - phiSector;
      ActsScalar phiMax = phiCenter + phiSector;
      // Fill the maps
      fillMap(rMap, {rMin, rMax});
      fillMap(zMap, {zMin, zMax});
      fillMap(phiMap, {phiMin, phiMax});
    }
  }

  for (auto [im, map] : enumerate(valueMaps)) {
    for (auto [key, value] : map) {
      boundaries[im].push_back(key);
    }
    std::sort(boundaries[im].begin(), boundaries[im].end());
  }

  ACTS_VERBOSE("- did yield " << boundaries[0u].size() << " boundaries in R.");
  ACTS_VERBOSE("- did yield " << boundaries[1u].size() << " boundaries in z.");
  ACTS_VERBOSE("- did yield " << boundaries[2u].size()
                              << " boundaries in phi.");

  return boundaries;
}
