// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/CylindricalDetectorHelper.hpp"

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/detail/PortalHelper.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

// Indexing of the portals follows the generation order of portals in the
// CylinderVolumeBounds and BevelledCylinderVolumeBounds (latter for wrapping)
//
// In short:
//
// 0: index of the low-z disc bound
// 1: index of the high-z disc bound
// 2: index of the outer cylinder bound
//
// If the volume doesn't extend inward to inner radius=0:
// 3: index of the inner cylinder bound
//
// Cylindrical sectors have up to 6 portals, enumerated as follows:
// 0: index of the low-z disc sector bound
// 1: index of the high-z disc sector bound
// 2: index of the outer cylinder sector bound
//
// If no inner cylinder sector bound is present:
// 3: index of the low-phi sector bound
// 4: index of the high-phi sector bound
// If an inner cylinder sector bound exists, it takes index 3 and the phi sector
// bounds are shifted by one index: 3: index of the inner cylinder sector bound
// 4: index of the low-phi sector bound
// 5: index of the high-phi sector bound
//

namespace {

/// @brief Helper function to create a disc portal replacement
///
/// @param transform The transform of the newly created portal
/// @param rBoundaries the vector of binning boundaries in r
/// @param phiBoundaries an eventual phi sector value
/// @param dir the  direction to be set
///
/// @return a new portal replacement object
Acts::Experimental::PortalReplacement createDiscReplacement(
    const Acts::Transform3& transform,
    const std::vector<Acts::ActsScalar>& rBoundaries,
    const std::vector<Acts::ActsScalar>& phiBoundaries, unsigned int index,
    Acts::Direction dir) {
  // Autodetector stitch value
  Acts::BinningValue stitchValue =
      phiBoundaries.size() == 2u ? Acts::binR : Acts::binPhi;
  // Estimate ranges
  auto [minR, maxR] = Acts::min_max(rBoundaries);
  auto [sectorPhi, avgPhi] = Acts::range_medium(phiBoundaries);

  // Transform and bounds
  auto bounds =
      std::make_unique<Acts::RadialBounds>(minR, maxR, 0.5 * sectorPhi, avgPhi);
  // A new surface on the negative side over the full range
  auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(
      transform, std::move(bounds));
  // Make a portal and indicate the new link direction
  const auto& stitchBoundaries =
      (stitchValue == Acts::binR) ? rBoundaries : phiBoundaries;
  return Acts::Experimental::PortalReplacement(
      Acts::Experimental::Portal::makeShared(surface), index, dir,
      stitchBoundaries, stitchValue);
}

/// @brief Helper function to create a cylinder portal replacement
///
/// @param transform The transform of the newly created portal
/// @param r is the radius of the cylinder
/// @param zBoundaries the vector of binning boundaries
/// @param phiBoundaries the vector of binning boundaries
/// @param index the index of this portal to be set
/// @param dir the navigation direction to be set
///
/// @return a new portal replacement object
Acts::Experimental::PortalReplacement createCylinderReplacement(
    const Acts::Transform3& transform, Acts::ActsScalar r,
    const std::vector<Acts::ActsScalar>& zBoundaries,
    const std::vector<Acts::ActsScalar>& phiBoundaries, unsigned int index,
    Acts::Direction dir) {
  // Autodetector stitch value
  Acts::BinningValue stitchValue =
      phiBoundaries.size() == 2u ? Acts::binZ : Acts::binPhi;
  auto [lengthZ, medZ] = Acts::range_medium(zBoundaries);
  auto [sectorPhi, avgPhi] = Acts::range_medium(phiBoundaries);

  // New bounds, with current length and sector values
  auto bounds = std::make_unique<Acts::CylinderBounds>(r, 0.5 * lengthZ,
                                                       0.5 * sectorPhi, avgPhi);
  // A new surface on the negative side over the full range
  auto surface = Acts::Surface::makeShared<Acts::CylinderSurface>(
      transform, std::move(bounds));

  // A make a portal and indicate the new link direction
  const auto& stitchBoundaries =
      (stitchValue == Acts::binZ) ? zBoundaries : phiBoundaries;
  return Acts::Experimental::PortalReplacement(
      Acts::Experimental::Portal::makeShared(surface), index, dir,
      stitchBoundaries, stitchValue);
}

/// @brief Helper function to create a disc portal replacement
/// @param gcxt the geometry context of this call
/// @param volumeCenter a reference center of the volume (combined)
/// @param refSurface the reference surface (old portal)
/// @param boundaries the vector of binning boundaries in r
/// @param binning the binning of the sector (inR, inZ)
/// @param index the index of this portal to be set
/// @param dir the navigation direction to be set
///
/// @return a new portal replacement object
Acts::Experimental::PortalReplacement createSectorReplacement(
    const Acts::GeometryContext& gctx, const Acts::Vector3& volumeCenter,
    const Acts::Surface& refSurface,
    const std::vector<Acts::ActsScalar>& boundaries, Acts::BinningValue binning,
    unsigned int index, Acts::Direction dir) {
  // Get a reference transform
  const auto& refTransform = refSurface.transform(gctx);
  auto refRotation = refTransform.rotation();
  // Bounds handling
  const auto& boundValues = refSurface.bounds().values();
  std::unique_ptr<Acts::PlanarBounds> bounds = nullptr;

  // Create a new transform
  Acts::Transform3 transform = Acts::Transform3::Identity();
  if (binning == Acts::binR) {
    // Range and center-r calculation
    auto [range, medium] = Acts::range_medium(boundaries);
    // New joint center:
    // - you start from the center of the volume, and then head the distence
    //   of medium-r along the y-axis of the former (and new) portal
    Acts::Vector3 pCenter = volumeCenter + medium * refRotation.col(1u);
    transform.pretranslate(pCenter);
    // Create the halflength
    Acts::ActsScalar halfX =
        0.5 * (boundValues[Acts::RectangleBounds::BoundValues::eMaxX] -
               boundValues[Acts::RectangleBounds::BoundValues::eMinX]);
    // New joint bounds
    bounds = std::make_unique<Acts::RectangleBounds>(halfX, 0.5 * range);
  } else if (binning == Acts::binZ) {
    // Range and medium z alculation
    auto [range, medium] = Acts::range_medium(boundaries);
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
    bounds = std::make_unique<Acts::RectangleBounds>(0.5 * range, halfY);
  }
  // Set the rotation
  transform.prerotate(refRotation);
  // A new surface on the negative side over the full range
  auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      transform, std::move(bounds));
  // A make a portal and indicate the new link direction
  Acts::Experimental::PortalReplacement pRep = {
      Acts::Experimental::Portal::makeShared(surface), index, dir, boundaries,
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
    const std::vector<Acts::Experimental::DetectorComponent::PortalContainer>&
        containers,
    const std::vector<unsigned int>& sides,
    const std::vector<unsigned int>& selectedOnly = {},
    Acts::Logging::Level logLevel = Acts::Logging::INFO) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("::stripSideVolumes", logLevel));

  // Thes are the stripped off outside volumes
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

  // Loop through the containers
  for (const auto& pc : containers) {
    // Loop through the selected sides and check if they are contained
    for (const auto& s : selectedSides) {
      auto cSide = pc.find(s);
      if (cSide != pc.end()) {
        auto p = cSide->second;
        auto& sVolumes = sideVolumes[s];
        auto aVolumes =
            Acts::Experimental::detail::PortalHelper::attachedDetectorVolumes(
                *p);
        sVolumes.insert(sVolumes.end(), aVolumes.begin(), aVolumes.end());
      }
    }
  }
  // return them
  return sideVolumes;
}

/// @brief Helper method to check alignment of the volumes, this method checks
/// if the z-axes are aligned and throws an exception if not, it assumes that
/// the detector volume content has been checked already
///
/// @param gctx the geometry context
/// @param volumes the input volumes to be checked
///
/// @note this is a strict matching that requires the rotation to be identical
///
/// @note throws exception if any of checks fails
void checkAlignment(
    const Acts::GeometryContext& gctx,
    const std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>&
        volumes) {
  // Take first transform as reference transform
  auto refRotation = volumes[0u]->transform(gctx).rotation();
  // Loop over rest and recursively test
  for (auto [iv, v] : Acts::enumerate(volumes)) {
    if (iv > 0) {
      auto curRotation = v->transform(gctx).rotation();
      if (not curRotation.isApprox(refRotation)) {
        std::string message = "CylindricalDetectorHelper: rotation of volume ";
        message += std::to_string(iv);
        message += std::string(" is not aligned with previous volume");
        throw std::invalid_argument(message.c_str());
      }
    }
  }
}

/// @brief Helper method to check the volumes in general and throw and excpetion if failes
///
/// @param gctx the geometry context
/// @param volumes the input volumes to be checked
///
/// @note throws exception if any of checks fails
void checkVolumes(
    const Acts::GeometryContext& gctx,
    const std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>&
        volumes) {
  // A minimal set of checks - exceptions are thrown
  // - not enough volumes given
  std::string message = "CylindricalDetectorHelper: ";
  if (volumes.size() < 2u) {
    message += std::string("not enough volume given (") +
               std::to_string(volumes.size());
    message += std::string(" ), when required >=2.");
    throw std::invalid_argument(message.c_str());
  }
  // - null pointer detected or non-cylindrical volume detected
  for (auto [iv, v] : Acts::enumerate(volumes)) {
    // Check for nullptr
    if (v == nullptr) {
      message += std::string("nullptr detector instead of volume " +
                             std::to_string(iv));
      throw std::invalid_argument(message.c_str());
    }
    // Check for cylindrical volume type
    if (v->volumeBounds().type() != Acts::VolumeBounds::BoundsType::eCylinder) {
      message +=
          std::string("non-cylindrical volume bounds detected for volume " +
                      std::to_string(iv));
      throw std::invalid_argument(message.c_str());
    }
  }
  // Check the alignment of the volumes
  checkAlignment(gctx, volumes);
}

/// @brief Helper method to check the volume bounds
///
/// @param gctx the geometry context
/// @param volumes the input volumes to be checked
/// @param refCur reference versus current check
///
/// @note throws exception if any of checks fails
void checkBounds(
    [[maybe_unused]] const Acts::GeometryContext& gctx,
    const std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>&
        volumes,
    const std::vector<std::array<unsigned int, 2u>>& refCur) {
  // Reference values
  auto refValues = volumes[0u]->volumeBounds().values();
  for (auto [iv, v] : Acts::enumerate(volumes)) {
    if (iv > 0u) {
      auto curValues = v->volumeBounds().values();
      for (auto [r, c] : refCur) {
        if (std::abs(refValues[r] - curValues[c]) >
            Acts::s_onSurfaceTolerance) {
          std::string message = "CylindricalDetectorHelper: '";
          message += volumes[iv - 1]->name();
          if (r != c) {
            message += "' does not attach to '";
          } else {
            message += "' does not match with '";
          }
          message += volumes[iv]->name();
          message += "\n";
          message += " - at bound values ";
          message += std::to_string(refValues[r]);
          message += " / " + std::to_string(curValues[c]);
          throw std::runtime_error(message.c_str());
        }
      }
      refValues = curValues;
    }
  }
}

}  // namespace

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CylindricalDetectorHelper::connectInR(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<Experimental::DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) {
  // Basic checks for eligability of the volumes
  checkVolumes(gctx, volumes);
  // Check for bounds values of volumes (i) and (i+1) succesively for
  // compatibility:
  // - check outer (1u) of volume (i) vs inner radius (0u) of volume (i+1)
  // - phi sector (3u) and average phi (4u) between volumes (i), (i+1)
  std::vector<std::array<unsigned int, 2u>> checks = {
      {1u, 0u}, {3u, 3u}, {4u, 4u}};
  // - And we check the half length if it is not a selective attachment
  if (selectedOnly.empty()) {
    checks.push_back({2u, 2u});
  }
  checkBounds(gctx, volumes, checks);

  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << volumes.size() << " detector volumes in R.");

  // Return proto container
  DetectorComponent::PortalContainer dShell;

  // Innermost volume boundaries
  std::vector<ActsScalar> rBoundaries = {};
  auto refValues = volumes[0u]->volumeBounds().values();

  // Reference boundary values
  rBoundaries.push_back(refValues[CylinderVolumeBounds::BoundValues::eMinR]);
  rBoundaries.push_back(refValues[CylinderVolumeBounds::BoundValues::eMaxR]);

  // Connect in R ? (2u is the index of the outer cylinder)
  bool connectR = selectedOnly.empty() or
                  std::find(selectedOnly.begin(), selectedOnly.end(), 2u) !=
                      selectedOnly.end();

  // Get phi sector and average phi
  ActsScalar phiSector =
      refValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector];
  ActsScalar avgPhi = refValues[CylinderVolumeBounds::BoundValues::eAveragePhi];

  // Fuse the cylinders - portals can be reused for this operation
  for (unsigned int iv = 1; iv < volumes.size(); ++iv) {
    refValues = volumes[iv]->volumeBounds().values();
    // Keep on collecting the outside maximum r for the overal r boundaries
    rBoundaries.push_back(refValues[CylinderVolumeBounds::BoundValues::eMaxR]);
    // Only connect if configured to do so
    if (connectR) {
      ACTS_VERBOSE("Connect volume '" << volumes[iv - 1]->name() << "' to "
                                      << volumes[iv]->name() << "'.");

      // When fusing volumes at a cylinder boundary, we *keep* one
      // portal and tranfer the portal link information from the other
      //
      // In this case the outer cylinder portal of the inner volume is kept,
      // the inner cylinder of the outer portal goes to waste
      auto& keepCylinder = volumes[iv - 1]->portalPtrs()[2u];
      auto& wasteCylinder = volumes[iv]->portalPtrs()[3u];
      keepCylinder->fuse(wasteCylinder);
      volumes[iv]->updatePortal(keepCylinder, 3u);
    }
  }

  // Proto container setting
  if (connectR) {
    // The number of portals indicate again if the inner is present or not,
    if (volumes[0u]->portals().size() == 4u or
        volumes[0u]->portals().size() == 6u) {
      dShell[3u] = volumes[0u]->portalPtrs()[3u];
    }
    dShell[2u] = volumes[volumes.size() - 1u]->portalPtrs()[2u];
  }

  // Check if sectors are present by the number of portals, check is done on the
  // outermost volume as this is required to have an inner cylinder, and hence
  // no offset needs to be respected
  bool sectorsPresent = volumes[volumes.size() - 1u]->portals().size() > 4u;

  // A portal replacement, it comprises of the portal, the index, the
  // direction, the binning and bins
  std::vector<PortalReplacement> pReplacements = {};

  // Disc assignments are forwad for negative disc, backward for positive
  std::vector<Acts::Direction> discDirs = {Acts::Direction::Forward,
                                           Acts::Direction::Backward};
  for (const auto [iu, idir] : enumerate(discDirs)) {
    if (selectedOnly.empty() or
        std::find(selectedOnly.begin(), selectedOnly.end(), iu) !=
            selectedOnly.end()) {
      const Surface& refSurface = volumes[0u]->portals()[iu]->surface();
      const Transform3& refTransform = refSurface.transform(gctx);
      pReplacements.push_back(createDiscReplacement(
          refTransform, rBoundaries, {avgPhi - phiSector, avgPhi + phiSector},
          iu, idir));
    }
  }

  // If volume sectors are present, these have to be respected
  if (sectorsPresent) {
    ACTS_VERBOSE("Sector planes are present, they need replacement.");
    // Sector assignments are forward backward
    std::vector<Acts::Direction> sectorDirs = {Acts::Direction::Forward,
                                               Acts::Direction::Backward};
    Acts::Vector3 vCenter = volumes[0u]->transform(gctx).translation();
    for (const auto [iu, idir] : enumerate(sectorDirs)) {
      // (iu + 4u) corresponds to the indices of the phi-low and phi-high sector
      // planes.
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
  } else {
    ACTS_VERBOSE(
        "No sector planes present, full 2 * M_PI cylindrical geometry.");
  }

  // Attach the new volume multi links
  PortalHelper::attachDetectorVolumeUpdators(gctx, volumes, pReplacements);

  // Exchange the portals of the volumes
  ACTS_VERBOSE("Portals of " << volumes.size() << " volumes need updating.");
  // Exchange the portals of the volumes
  for (auto& iv : volumes) {
    ACTS_VERBOSE("- update portals of volume '" << iv->name() << "'.");
    for (auto& [p, i, dir, boundaries, binning] : pReplacements) {
      // Fill the map
      dShell[i] = p;

      // Potential offset for tube vs/ cylinder
      // - if the volume doesn't have an inner portal, indices need to
      //   be shifted by -1 to update the correct index, that's the case for
      //   size 3 and 5 for portals
      size_t nPortals = iv->portals().size();
      bool innerPresent = (nPortals == 3u or nPortals == 5u);
      int iOffset = (innerPresent and i > 2u) ? -1 : 0;
      ACTS_VERBOSE("-- update portal with index "
                   << i + iOffset << " (including offset " << iOffset << ")");
      iv->updatePortal(p, static_cast<unsigned int>(i + iOffset));
    }
  }
  // Done.
  return dShell;
}

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CylindricalDetectorHelper::connectInZ(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>& volumes,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) {
  // Basic checks for eligability of the volumes
  checkVolumes(gctx, volumes);
  // Check for bounds compatibility
  // We check phi sector (3u) and average phi (4u)
  std::vector<std::array<unsigned int, 2u>> checks = {{3u, 3u}, {4u, 4u}};
  // And we check the inner radius [0u], outer radius[1u] if its' not a
  // selective attachment
  if (selectedOnly.empty()) {
    checks.push_back({0u, 0u});
    checks.push_back({1u, 1u});
  }
  checkBounds(gctx, volumes, checks);

  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << volumes.size() << " detector volumes in Z.");

  // Return proto container
  DetectorComponent::PortalContainer dShell;

  // Connect in Z ?
  // - 1u corresponds to the index of the high-z disc portal for the reference
  // volume.
  bool connectZ = selectedOnly.empty() or
                  std::find(selectedOnly.begin(), selectedOnly.end(), 1u) !=
                      selectedOnly.end();
  // Reference z axis
  const auto rotation = volumes[0u]->transform(gctx).rotation();

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
      ACTS_VERBOSE("Connect volume '" << volumes[iv - 1]->name() << "' to "
                                      << volumes[iv]->name() << "'.");
      // When fusing, one portal survives (keep) and gets the
      // portal linking from the waste tranfered
      //
      // In this case we keep the disc at positive z of the volume
      // at lower relative z, and trash the disc at negative z of the
      // following volume
      auto& keepDisc = volumes[iv - 1]->portalPtrs()[1u];
      auto& wasteDisc = volumes[iv]->portalPtrs()[0u];
      // Throw an exception if the discs are not at the same position
      Vector3 keepPosition = keepDisc->surface().center(gctx);
      Vector3 wastePosition = wasteDisc->surface().center(gctx);
      if (not keepPosition.isApprox(wastePosition)) {
        std::string message = "CylindricalDetectorHelper: '";
        message += volumes[iv - 1]->name();
        message += "' does not attach to '";
        message += volumes[iv]->name();
        message += "\n";
        message += " - along z with values ";
        message += Acts::toString(keepPosition);
        message += " / " + Acts::toString(wastePosition);
        throw std::runtime_error(message.c_str());
      }
      keepDisc->fuse(wasteDisc);
      volumes[iv]->updatePortal(keepDisc, 0u);
    }
  }

  // Register what we have from the container
  if (connectZ) {
    dShell[0u] = volumes[0u]->portalPtrs()[0u];
    dShell[1u] = volumes[volumes.size() - 1u]->portalPtrs()[1u];
  }

  // Centre of gravity
  Vector3 combinedCenter =
      0.5 * (zBoundaries3D[zBoundaries3D.size() - 1u] + zBoundaries3D[0u]);

  ACTS_VERBOSE("New combined center calculated at "
               << toString(combinedCenter));

  // Evaluate the series of z boundaries
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

  // A portal replacement, it comprises of the portal, the index, the
  // direction, the binning and bins
  std::vector<PortalReplacement> pReplacements = {};

  // Disc assignments are forwad for negative disc, backward for positive
  std::vector<Acts::Direction> cylinderDirs = {Acts::Direction::Backward};
  // Cylinder radii
  std::vector<Acts::ActsScalar> cylinderR = {maxR};
  if (innerPresent) {
    ACTS_VERBOSE("Inner surface present, tube geometry detected.");
    cylinderDirs.push_back(Direction::Forward);
    cylinderR.push_back(minR);
  } else {
    ACTS_VERBOSE("No inner surface present, solid cylinder geometry detected.");
  }
  // Tube/cylinder offset
  unsigned int iSecOffset = innerPresent ? 4u : 3u;
  // Prepare the cylinder replacements
  for (const auto [iu, idir] : enumerate(cylinderDirs)) {
    if (selectedOnly.empty() or
        std::find(selectedOnly.begin(), selectedOnly.end(), iu + 2u) !=
            selectedOnly.end()) {
      pReplacements.push_back(createCylinderReplacement(
          combinedTransform, cylinderR[iu], zBoundaries,
          {avgPhi - phiSector, avgPhi + phiSector}, iu + 2u, idir));
    }
  }

  // Prepare the sector side replacements
  if (sectorsPresent) {
    ACTS_VERBOSE("Sector planes are present, they need replacement.");
    // Sector assignmenta are forward backward
    std::vector<Acts::Direction> sectorDirs = {Acts::Direction::Forward,
                                               Acts::Direction::Backward};
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
  } else {
    ACTS_VERBOSE(
        "No sector planes present, full 2 * M_PI cylindrical geometry.");
  }

  // Attach the new volume multi links
  PortalHelper::attachDetectorVolumeUpdators(gctx, volumes, pReplacements);

  // Exchange the portals of the volumes
  ACTS_VERBOSE("Portals of " << volumes.size() << " volumes need updating.");
  for (auto& iv : volumes) {
    ACTS_VERBOSE("- update portals of volume '" << iv->name() << "'.");
    for (auto& [p, i, dir, boundaries, binning] : pReplacements) {
      // Potential offset for tube vs/ cylinder
      // if the volume doesn't have an inner portal, indices need to be shifted
      // by -1 to update the correct index.
      int iOffset = (i > 2u and not innerPresent) ? -1 : 0;
      ACTS_VERBOSE("-- update portal with index " << i);
      iv->updatePortal(p, static_cast<unsigned int>(i + iOffset));
      // Fill the map
      dShell[i] = p;
    }
  }
  // Done.
  return dShell;
}

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CylindricalDetectorHelper::connectInPhi(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>& volumes,
    const std::vector<unsigned int>& /*selectedOnly*/,
    Acts::Logging::Level logLevel) {
  // Basic checks for eligability of the volumes
  checkVolumes(gctx, volumes);

  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << volumes.size() << " detector volumes in phi.");

  // Return proto container
  DetectorComponent::PortalContainer dShell;

  // Check if inner cylinder and sectors are present by the number of portals
  size_t nPortals = volumes[volumes.size() - 1u]->portals().size();
  bool innerPresent = (nPortals != 3u and nPortals != 5u);

  Transform3 refTransform = volumes[0u]->transform(gctx);

  // Sector offset
  unsigned int iSecOffset = innerPresent ? 4u : 3u;
  std::vector<ActsScalar> phiBoundaries = {};
  auto refValues = volumes[0u]->volumeBounds().values();
  phiBoundaries.push_back(
      refValues[CylinderVolumeBounds::BoundValues::eAveragePhi] -
      refValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector]);
  phiBoundaries.push_back(
      refValues[CylinderVolumeBounds::BoundValues::eAveragePhi] +
      refValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector]);
  // Fuse the sectors - portals can be reused for this operation
  for (unsigned int iv = 1; iv < volumes.size(); ++iv) {
    ACTS_VERBOSE("Connect volume '" << volumes[iv - 1]->name() << "' to "
                                    << volumes[iv]->name() << "'.");

    // Fuse and swap
    auto& keepSector = volumes[iv - 1]->portalPtrs()[iSecOffset + 1u];
    auto& wasteSector = volumes[iv]->portalPtrs()[iSecOffset];
    keepSector->fuse(wasteSector);
    volumes[iv]->updatePortal(keepSector, iSecOffset);
    // The current values
    auto curValues = volumes[iv]->volumeBounds().values();
    // Bail out if they do not match
    ActsScalar lowPhi =
        curValues[CylinderVolumeBounds::BoundValues::eAveragePhi] -
        curValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector];
    ActsScalar highPhi =
        curValues[CylinderVolumeBounds::BoundValues::eAveragePhi] +
        curValues[CylinderVolumeBounds::BoundValues::eHalfPhiSector];
    // Check phi attachment
    if (std::abs(phiBoundaries[phiBoundaries.size() - 1u] - lowPhi) >
        Acts::s_onSurfaceTolerance) {
      std::string message = "CylindricalDetectorHelper: '";
      message += volumes[iv - 1]->name();
      message += "' does not attach to '";
      message += volumes[iv]->name();
      message += "\n";
      message += " - within phi sectors ";
      message += std::to_string(lowPhi);
      message +=
          " / " + std::to_string(phiBoundaries[phiBoundaries.size() - 1u]);
      throw std::runtime_error(message.c_str());
    }
    // Check radial and longitudinal compatibility - @TODO
    phiBoundaries.push_back(highPhi);
    // Recursive setting of the values
    refValues = curValues;
  }

  // A portal replacement, it comprises of the portal, the index, the
  // direction, the binning and bins
  std::vector<PortalReplacement> pReplacements = {};
  // Negative disc
  pReplacements.push_back(createDiscReplacement(
      refTransform,
      {refValues[CylinderVolumeBounds::BoundValues::eMinR],
       refValues[CylinderVolumeBounds::BoundValues::eMaxR]},
      phiBoundaries, 0u, Acts::Direction::Forward));

  // Positive disc
  pReplacements.push_back(createDiscReplacement(
      refTransform,
      {refValues[CylinderVolumeBounds::BoundValues::eMinR],
       refValues[CylinderVolumeBounds::BoundValues::eMaxR]},
      phiBoundaries, 1u, Acts::Direction::Backward));

  // Outside cylinder
  pReplacements.push_back(createCylinderReplacement(
      refTransform, refValues[CylinderVolumeBounds::BoundValues::eMaxR],
      {-refValues[CylinderVolumeBounds::BoundValues::eHalfLengthZ],
       refValues[CylinderVolumeBounds::BoundValues::eHalfLengthZ]},
      phiBoundaries, 2u, Acts::Direction::Backward));

  // If the volume has a different inner radius than 0, it MUST have
  // an inner cylinder
  if (refValues[CylinderVolumeBounds::BoundValues::eMinR] > 0.) {
    // Inner cylinder
    pReplacements.push_back(createCylinderReplacement(
        refTransform, refValues[CylinderVolumeBounds::BoundValues::eMinR],
        {-refValues[CylinderVolumeBounds::BoundValues::eHalfLengthZ],
         refValues[CylinderVolumeBounds::BoundValues::eHalfLengthZ]},
        phiBoundaries, 3u, Acts::Direction::Forward));
  }

  // Attach the new volume multi links
  PortalHelper::attachDetectorVolumeUpdators(gctx, volumes, pReplacements);
  // Exchange the portals of the volumes
  ACTS_VERBOSE("Portals of " << volumes.size() << " volumes need updating.");
  for (auto& iv : volumes) {
    ACTS_VERBOSE("- update portals of volume '" << iv->name() << "'.");
    for (auto& [p, i, dir, boundaries, binning] : pReplacements) {
      ACTS_VERBOSE("-- update portal with index " << i);
      iv->updatePortal(p, static_cast<unsigned int>(i));
      // Fill the map
      dShell[i] = p;
    }
  }

  // Done.
  return dShell;
}

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CylindricalDetectorHelper::wrapInZR(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>>& volumes,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Wrapping volumes in Z-R.");

  // Minimal set of checks
  if (volumes.size() != 2u) {
    throw std::invalid_argument(
        "Wrapping the detector volume requires exactly 2 volumes.");
  }

  // Return the new container
  DetectorComponent::PortalContainer dShell;

  // Keep the outer shells
  dShell[0u] = volumes[1u]->portalPtrs()[0u];
  dShell[1u] = volumes[1u]->portalPtrs()[1u];
  dShell[2u] = volumes[1u]->portalPtrs()[2u];

  // Fuse outer cover of first with inner cylinder of wrapping volume
  auto& keepCover = volumes[0u]->portalPtrs()[2u];
  auto& wasteCover = volumes[1u]->portalPtrs()[3u];
  keepCover->fuse(wasteCover);
  volumes[1u]->updatePortal(keepCover, 3u);

  // Stitch sides - negative
  auto& keepDiscN = volumes[1u]->portalPtrs()[4u];
  auto& wasteDiscN = volumes[0u]->portalPtrs()[0u];
  keepDiscN->fuse(wasteDiscN);
  volumes[0u]->updatePortal(keepDiscN, 0u);

  // Stich sides - positive
  auto& keepDiscP = volumes[0u]->portalPtrs()[1u];
  auto& wasteDiscP = volumes[1u]->portalPtrs()[5u];
  keepDiscP->fuse(wasteDiscP);
  volumes[1u]->updatePortal(keepDiscP, 5u);

  // If needed, insert new cylinder
  if (volumes[0u]->portalPtrs().size() == 4u and
      volumes[1u]->portalPtrs().size() == 8u) {
    // We need a new cylinder spanning over the entire inner tube
    ActsScalar hlZ =
        volumes[0u]
            ->volumeBounds()
            .values()[Acts::CylinderVolumeBounds::BoundValues::eHalfLengthZ];
    ActsScalar HlZ =
        volumes[1u]->volumeBounds().values()
            [Acts::CutoutCylinderVolumeBounds::BoundValues::eHalfLengthZ];
    ActsScalar innerR =
        volumes[0u]
            ->volumeBounds()
            .values()[Acts::CylinderVolumeBounds::BoundValues::eMinR];
    // Create the inner replacement
    std::vector<PortalReplacement> pReplacements;
    pReplacements.push_back(createCylinderReplacement(
        volumes[0u]->transform(gctx), innerR, {-HlZ, -hlZ, hlZ, HlZ},
        {-M_PI, M_PI}, 3u, Direction::Forward));
    std::vector<std::shared_ptr<DetectorVolume>> zVolumes = {
        volumes[1u], volumes[0u], volumes[1u]};
    // Attach the new volume multi links
    PortalHelper::attachDetectorVolumeUpdators(gctx, zVolumes, pReplacements);
    auto& [p, i, dir, boundaries, binning] = pReplacements[0u];
    // Update the portals
    volumes[1u]->updatePortal(p, 6u);
    volumes[0u]->updatePortal(p, 3u);
    volumes[1u]->updatePortal(p, 7u);
    // Inner skin
    dShell[3u] = p;
  }
  // Done.
  return dShell;
}

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CylindricalDetectorHelper::connectInR(
    const GeometryContext& gctx,
    const std::vector<DetectorComponent::PortalContainer>& containers,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) noexcept(false) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << containers.size() << " proto containers in R.");

  // Return the new container
  DetectorComponent::PortalContainer dShell;

  // Fuse the cylinders - portals can be reused for this operation
  for (unsigned int ic = 1; ic < containers.size(); ++ic) {
    auto& formerContainer = containers[ic - 1];
    auto& currentContainer = containers[ic];
    // Check and throw exception
    if (formerContainer.find(2u) == formerContainer.end()) {
      throw std::invalid_argument(
          "CylindricalDetectorHelper: proto container has no outer cover, "
          "can "
          "not be connected in R");
    }
    if (currentContainer.find(3u) == currentContainer.end()) {
      throw std::invalid_argument(
          "CylindricalDetectorHelper: proto container has no inner cover, "
          "can "
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
    dShell[3u] = containers[0u].find(3u)->second;
  }
  dShell[2u] = containers[containers.size() - 1u].find(2u)->second;

  auto sideVolumes =
      stripSideVolumes(containers, {0u, 1u, 4u, 5u}, selectedOnly, logLevel);

  for (auto [s, volumes] : sideVolumes) {
    auto pR = connectInR(gctx, volumes, {s});
    if (pR.find(s) != pR.end()) {
      dShell[s] = pR.find(s)->second;
    }
  }

  // Done.
  return dShell;
}

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CylindricalDetectorHelper::connectInZ(
    const GeometryContext& gctx,
    const std::vector<DetectorComponent::PortalContainer>& containers,
    const std::vector<unsigned int>& selectedOnly,
    Acts::Logging::Level logLevel) noexcept(false) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Connect " << containers.size() << " proto containers in Z.");

  // Return the new container
  DetectorComponent::PortalContainer dShell;

  for (unsigned int ic = 1; ic < containers.size(); ++ic) {
    auto& formerContainer = containers[ic - 1];
    auto& currentContainer = containers[ic];
    // Check and throw exception
    if (formerContainer.find(1u) == formerContainer.end()) {
      throw std::invalid_argument(
          "CylindricalDetectorHelper: proto container has no negative disc, "
          "can not be connected in Z");
    }
    if (currentContainer.find(0u) == currentContainer.end()) {
      throw std::invalid_argument(
          "CylindricalDetectorHelper: proto container has no positive disc, "
          "can not be connected in Z");
    }
    std::shared_ptr<Portal> keepDisc = formerContainer.find(1u)->second;
    std::shared_ptr<Portal> wasteDisc = currentContainer.find(0u)->second;
    keepDisc->fuse(wasteDisc);
    for (auto& av : wasteDisc->attachedDetectorVolumes()[1u]) {
      ACTS_VERBOSE("Update portal of detector volume '" << av->name() << "'.");
      av->updatePortal(keepDisc, 0u);
    }
  }

  // Proto container refurbishment
  dShell[0u] = containers[0u].find(0u)->second;
  dShell[1u] = containers[containers.size() - 1u].find(1u)->second;

  // Check if this is a tube or a cylinder container (check done on 1st)
  std::vector<unsigned int> nominalSides = {2u, 4u, 5u};
  if (containers[0u].find(3u) != containers[0u].end()) {
    nominalSides.push_back(3u);
  }

  // Strip the side volumes
  auto sideVolumes =
      stripSideVolumes(containers, nominalSides, selectedOnly, logLevel);

  ACTS_VERBOSE("There remain " << sideVolumes.size()
                               << " side volume packs to be connected");
  for (auto [s, volumes] : sideVolumes) {
    ACTS_VERBOSE(" - connect " << volumes.size() << " at selected side " << s);
    auto pR = connectInZ(gctx, volumes, {s}, logLevel);
    if (pR.find(s) != pR.end()) {
      dShell[s] = pR.find(s)->second;
    }
  }

  // Done.
  return dShell;
}

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CylindricalDetectorHelper::connectInPhi(
    [[maybe_unused]] const GeometryContext& gctx,
    [[maybe_unused]] const std::vector<DetectorComponent::PortalContainer>&
        containers,
    [[maybe_unused]] const std::vector<unsigned int>& selectedOnly,
    [[maybe_unused]] Acts::Logging::Level logLevel) noexcept(false) {
  throw std::invalid_argument(
      "CylindricalDetectorHelper: container connection in phi not implemented "
      "yet.");
  DetectorComponent::PortalContainer dShell;
  // Done.
  return dShell;
}

Acts::Experimental::DetectorComponent::PortalContainer
Acts::Experimental::detail::CylindricalDetectorHelper::wrapInZR(
    [[maybe_unused]] const GeometryContext& gctx,
    const std::vector<DetectorComponent::PortalContainer>& containers,
    Acts::Logging::Level logLevel) {
  if (containers.size() != 2u) {
    throw std::invalid_argument(
        "CylindricalDetectorHelper: wrapping must take exaclty two "
        "containers.");
  }

  // The inner one is a container
  auto innerContainer = containers.front();
  // The outer one is a single volume represented as a container
  auto outerContainer = containers.back();
  std::shared_ptr<DetectorVolume> wrappingVolume = nullptr;
  for (auto [key, value] : outerContainer) {
    auto attachedVolumes = value->attachedDetectorVolumes();
    for (const auto& ava : attachedVolumes) {
      for (const auto& av : ava) {
        if (wrappingVolume == nullptr and av != nullptr) {
          wrappingVolume = av;
        } else if (wrappingVolume != nullptr and av != wrappingVolume) {
          throw std::invalid_argument(
              "CylindricalDetectorHelper: wrapping container must represent a "
              "single volume.");
        }
      }
    }
  }
  if (wrappingVolume == nullptr) {
    throw std::invalid_argument(
        "CylindricalDetectorHelper: wrapping volume could not be "
        "determined.");
  }

  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("CylindricalDetectorHelper", logLevel));

  ACTS_DEBUG("Wrapping a container with volume `" << wrappingVolume->name()
                                                  << "'.");
  // Return the new container
  DetectorComponent::PortalContainer dShell;

  // Keep the outer shells of the proto container
  dShell[0u] = wrappingVolume->portalPtrs()[0u];
  dShell[1u] = wrappingVolume->portalPtrs()[1u];
  dShell[2u] = wrappingVolume->portalPtrs()[2u];

  // Fuse outer cover of first with inner cylinder of wrapping volume
  auto& keepCover = innerContainer[2u];
  auto& wasteCover = wrappingVolume->portalPtrs()[3u];
  keepCover->fuse(wasteCover);
  wrappingVolume->updatePortal(keepCover, 3u);

  // Stitch sides - negative
  auto& keepDiscN = innerContainer[0u];
  auto& wasteDiscN = wrappingVolume->portalPtrs()[4u];
  keepDiscN->fuse(wasteDiscN);
  wrappingVolume->updatePortal(keepDiscN, 4u);

  // Stich sides - positive
  auto& keepDiscP = innerContainer[1u];
  auto& wasteDiscP = wrappingVolume->portalPtrs()[5u];
  keepDiscP->fuse(wasteDiscP);
  wrappingVolume->updatePortal(keepDiscP, 5u);

  // If inner stitching is necessary
  if (innerContainer.size() == 4u and
      wrappingVolume->portalPtrs().size() == 8u) {
    // Inner Container portal
    auto& centralSegment = innerContainer[3u];
    auto centralValues = centralSegment->surface().bounds().values();
    ActsScalar centralHalfLengthZ =
        centralValues[CylinderBounds::BoundValues::eHalfLengthZ];
    // The two segments
    auto& nSegment = wrappingVolume->portalPtrs()[6u];
    auto nValues = nSegment->surface().bounds().values();
    ActsScalar nHalfLengthZ =
        nValues[CylinderBounds::BoundValues::eHalfLengthZ];
    auto& pSegment = wrappingVolume->portalPtrs()[7u];
    auto pValues = pSegment->surface().bounds().values();
    ActsScalar pHalfLengthZ =
        pValues[CylinderBounds::BoundValues::eHalfLengthZ];

    auto sideVolumes = stripSideVolumes({innerContainer}, {3u}, {3u}, logLevel);

    // First the left volume sector
    std::vector<std::shared_ptr<DetectorVolume>> innerVolumes = {
        wrappingVolume->getSharedPtr()};

    std::vector<ActsScalar> zBoundaries = {
        -centralHalfLengthZ - 2 * nHalfLengthZ, centralHalfLengthZ};
    // Loop over side volume and register the z boundaries
    for (auto& svs : sideVolumes) {
      for (auto& v : svs.second) {
        ActsScalar hlZ = v->volumeBounds().values()[2u];
        zBoundaries.push_back(zBoundaries.back() + 2 * hlZ);
        innerVolumes.push_back(v);
      }
    }
    // Last the right volume sector
    zBoundaries.push_back(zBoundaries.back() + 2 * pHalfLengthZ);
    innerVolumes.push_back(wrappingVolume);
  }

  // Done.
  return dShell;
}

std::array<std::vector<Acts::ActsScalar>, 3u>
Acts::Experimental::detail::CylindricalDetectorHelper::rzphiBoundaries(
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
