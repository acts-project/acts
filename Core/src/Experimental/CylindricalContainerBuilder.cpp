// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/CylindricalContainerBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Enumerate.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/PortalLinks.hpp"
#include "Acts/Experimental/VolumeLinks.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/detail/Axis.hpp"

#include <memory>
#include <stdexcept>

std::shared_ptr<Acts::DetectorVolume> Acts::CylindricalContainerBuilderR::build(
    const GeometryContext& gctx, const VolumeBuilder& volumeBuilder,
    const std::string& name) {
  auto containerVolumes = volumeBuilder(gctx);

  // Consistency check: volumes exist
  if (containerVolumes.empty()) {
    return nullptr;
  }

  // Consistency check: positioning, alignment, bounds
  //
  // & remember the ranges and borders

  // The first volume
  DetectorVolume* fVolume = containerVolumes.begin()->get();
  const Transform3& fTransform = fVolume->transform();
  std::vector<ActsScalar> fValues = fVolume->volumeBounds().values();

  DetectorVolume* lVolume = fVolume;

  // The values that need to correspond
  std::vector<CylinderVolumeBounds::BoundValues> checkValues = {
      CylinderVolumeBounds::eHalfLengthZ, CylinderVolumeBounds::eHalfPhiSector,
      CylinderVolumeBounds::eAveragePhi};

  std::vector<ActsScalar> rBoundaries = {};

  std::vector<PortalLink> nPortalLinks = {};
  std::vector<PortalLink> pPortalLinks = {};

  // Loop through, check & remember, exchange boundaries & collect links
  for (auto [i, v] : enumerate(containerVolumes)) {
    // The current volume
    DetectorVolume* cVolume = v.get();
    const Transform3& cTransform = cVolume->transform();
    std::vector<ActsScalar> cValues = cVolume->volumeBounds().values();

    auto cPortals = cVolume->portals();
    nPortalLinks.push_back(cPortals[0]->portalLink(forward));
    pPortalLinks.push_back(cPortals[1]->portalLink(backward));

    // The first one is further excluded from cross-comparison
    if (i == 0) {
      rBoundaries.push_back(cValues[CylinderVolumeBounds::eMinR]);
      rBoundaries.push_back(cValues[CylinderVolumeBounds::eMaxR]);
      continue;
    }
    rBoundaries.push_back(cValues[CylinderVolumeBounds::eMaxR]);

    // The rotation needs to be aligned
    if (not fTransform.rotation().isApprox(cTransform.rotation())) {
      throw std::invalid_argument(
          "\n *** Cylindrical container (in R): rotations are not aligned.");
    }

    // The half length, phi sector, etc, need to correspond
    for (auto cv : checkValues) {
      if (std::abs(fValues[cv] - cValues[cv]) >
          std::numeric_limits<ActsScalar>::epsilon()) {
        std::string sException(
            "\n *** Cylindrical container (in R): bounds are not "
            "consistent.\n");
        sException += std::string(" *** Failed on comparing ");
        sException += std::to_string(fValues[cv]);
        sException += " with ";
        sException += std::to_string(cValues[cv]);
        sException += " for parameter: ";
        sException += std::to_string(cv);
        throw std::invalid_argument(sException.c_str());
      }
    }

    // Update the portal
    auto cPortalPtrs = cVolume->portalPtrs();
    lVolume->updatePortalPtr(cPortalPtrs[3u], 2u, true, backward);
  }

  // Create the discs on left and right
  Vector3 center = fTransform.translation();
  Vector3 zAxis = fTransform.rotation().col(2);

  // Min R / Max R / half Phi and average Phi
  ActsScalar minR = rBoundaries[0];
  ActsScalar maxR = rBoundaries[rBoundaries.size() - 1];
  ActsScalar halfLengthZ = fValues[CylinderVolumeBounds::eHalfLengthZ];
  ActsScalar halfPhiSector = fValues[CylinderVolumeBounds::eHalfPhiSector];
  ActsScalar averagePhi = fValues[CylinderVolumeBounds::eAveragePhi];

  Vector3 nDiscCenter =
      center - fValues[CylinderVolumeBounds::eHalfLengthZ] * zAxis;
  Vector3 pDiscCenter =
      center + fValues[CylinderVolumeBounds::eHalfLengthZ] * zAxis;

  Transform3 nDiscTransform = Transform3::Identity();
  nDiscTransform.prerotate(fTransform.rotation());
  nDiscTransform.pretranslate(nDiscCenter);
  auto nDiscBounds =
      std::make_unique<RadialBounds>(minR, maxR, halfPhiSector, averagePhi);

  auto nDisc =
      Surface::makeShared<DiscSurface>(nDiscTransform, std::move(nDiscBounds));
  auto nPortal = std::make_shared<Portal>(std::move(nDisc));

  Transform3 pDiscTransform = Transform3::Identity();
  pDiscTransform.prerotate(fTransform.rotation());
  pDiscTransform.pretranslate(pDiscCenter);
  auto pDiscBounds =
      std::make_unique<RadialBounds>(minR, maxR, halfPhiSector, averagePhi);

  auto pDisc =
      Surface::makeShared<DiscSurface>(pDiscTransform, std::move(pDiscBounds));
  auto pPortal = std::make_shared<Portal>(std::move(pDisc));

  // A r-binned varialbe link: make three as we move
  VariableVolumeLink vrLink(detail::VariableAxis(rBoundaries), binR);
  VariableVolumeLink nrLink(detail::VariableAxis(rBoundaries), binR);
  VariableVolumeLink prLink(detail::VariableAxis(rBoundaries), binR);
  MultiplePortalLink nPortalMultiLink{nPortalLinks, std::move(nrLink)};
  MultiplePortalLink pPortalMultiLink{pPortalLinks, std::move(prLink)};

  nPortal->updatePortalLink(nPortalMultiLink, forward);
  pPortal->updatePortalLink(pPortalMultiLink, backward);

  // Update to the new portal surface
  for (auto [i, cv] : enumerate(containerVolumes)) {
    // The side portals
    cv->updatePortalPtr(nPortal, 0u, false);
    cv->updatePortalPtr(pPortal, 1u, false);
  }

  // Create new container bounds
  auto containerBounds = std::make_unique<CylinderVolumeBounds>(
      minR, maxR, halfLengthZ, halfPhiSector, averagePhi);

  return DetectorVolume::makeShared(fTransform, std::move(containerBounds),
                                    std::move(containerVolumes),
                                    std::move(vrLink), true, binR, name);
}

std::shared_ptr<Acts::DetectorVolume> Acts::CylindricalContainerBuilderZ::build(
    const GeometryContext& gctx, const VolumeBuilder& volumeBuilder,
    const std::string& name) {
  auto containerVolumes = volumeBuilder(gctx);
  // Consistency check: volumes exist
  if (containerVolumes.empty()) {
    return nullptr;
  }

  // Consistency check: positioning, alignment, bounds
  //
  // & remember the ranges and borders

  // The first volume
  DetectorVolume* fVolume = containerVolumes.begin()->get();
  const Transform3& fTransform = fVolume->transform();
  Vector3 fRef = fTransform.translation();

  Vector3 fAxisZ = fTransform.rotation().col(2);
  std::vector<ActsScalar> fValues = fVolume->volumeBounds().values();

  // The values that need to correspond
  std::vector<CylinderVolumeBounds::BoundValues> checkValues = {
      CylinderVolumeBounds::eMinR, CylinderVolumeBounds::eMaxR,
      CylinderVolumeBounds::eHalfPhiSector, CylinderVolumeBounds::eAveragePhi};

  ActsScalar zLengthTotal = 0.;
  std::vector<ActsScalar> zLengths = {};

  std::vector<PortalLink> oPortalLinks = {};
  std::vector<PortalLink> iPortalLinks = {};

  // Loop through, check and remember
  for (auto [i, v] : enumerate(containerVolumes)) {
    // The current volume
    DetectorVolume* cVolume = v.get();
    const Transform3& cTransform = cVolume->transform();
    auto cPortals = cVolume->portals();
    // Outer cover always exists
    oPortalLinks.push_back(cPortals[2]->portalLink(backward));
    // Inner portal might exist, if there's exactly 4 or 6 portals
    if (cPortals.size() == 4u or cPortals.size() == 6u) {
      iPortalLinks.push_back(cPortals[3]->portalLink(forward));
    }
    // TODO add sectoral portals

    std::vector<ActsScalar> cValues = cVolume->volumeBounds().values();
    zLengths.push_back(2 * cValues[CylinderVolumeBounds::eHalfLengthZ]);
    zLengthTotal += 2 * cValues[CylinderVolumeBounds::eHalfLengthZ];

    // The first one is further excluded from cross-comparison
    if (i == 0) {
      fRef = fRef - cValues[CylinderVolumeBounds::eHalfLengthZ] * fAxisZ;
      continue;
    }

    // The rotation needs to be aligned
    if (not fTransform.rotation().isApprox(cTransform.rotation())) {
      throw std::invalid_argument(
          "\n *** Cylindrical container (in Z): rotations are not aligned.");
    }

    // The inner/outer radii and the phiSector need to correspond
    for (auto cv : checkValues) {
      if (std::abs(fValues[cv] - cValues[cv]) >
          std::numeric_limits<ActsScalar>::epsilon()) {
        std::string sException(
            "\n *** Cylindrical container (in Z): bounds are not "
            "consistent.\n");
        sException += std::string(" *** Failed on comparing ");
        sException += std::to_string(fValues[cv]);
        sException += " with ";
        sException += std::to_string(cValues[cv]);
        sException += " for parameter: ";
        sException += std::to_string(cv);
        throw std::invalid_argument(sException.c_str());
      }
    }
  }  // end of loop

  // Create the z boundaries
  std::vector<ActsScalar> zBoundaries = {-0.5 * zLengthTotal};
  for (auto [i, zl] : enumerate(zLengths)) {
    zBoundaries.push_back(zBoundaries[i] + zl);
  }

  // New halflength
  fValues[CylinderVolumeBounds::eHalfLengthZ] = 0.5 * zLengthTotal;
  std::array<ActsScalar, CylinderVolumeBounds::eSize> containerValues;
  for (auto [i, v] : enumerate(fValues)) {
    containerValues[i] = v;
  }

  // The new container position
  Vector3 containerPosition = fRef + 0.5 * zLengthTotal * fAxisZ;
  auto containerRotation = fTransform.rotation();

  Transform3 containerTransform = Transform3::Identity();
  containerTransform.prerotate(containerRotation);
  containerTransform.pretranslate(containerPosition);

  // The outer, new, shared cylinder portal (always to be done)
  std::array<ActsScalar, CylinderBounds::eSize> outerValues = {
      fValues[CylinderVolumeBounds::eMaxR],
      fValues[CylinderVolumeBounds::eHalfLengthZ],
      fValues[CylinderVolumeBounds::eHalfPhiSector],
      fValues[CylinderVolumeBounds::eAveragePhi]};
  auto outerBounds = std::make_shared<const CylinderBounds>(outerValues);
  auto outerSurface =
      Surface::makeShared<CylinderSurface>(containerTransform, outerBounds);
  std::shared_ptr<Portal> outerPortal = std::make_shared<Portal>(outerSurface);
  // The inner one (is optional)
  std::shared_ptr<Portal> innerPortal = nullptr;
  if (fVolume->portals().size() == 4u) {
    std::array<ActsScalar, CylinderBounds::eSize> innerValues = {
        fValues[CylinderVolumeBounds::eMinR],
        fValues[CylinderVolumeBounds::eHalfLengthZ],
        fValues[CylinderVolumeBounds::eHalfPhiSector],
        fValues[CylinderVolumeBounds::eAveragePhi]};
    auto innerBounds = std::make_shared<const CylinderBounds>(innerValues);
    auto innerSurface =
        Surface::makeShared<CylinderSurface>(containerTransform, innerBounds);
    innerPortal = std::make_shared<Portal>(innerSurface);
  }

  // Set new portal surfaces
  for (auto v : containerVolumes) {
    v->updatePortalPtr(outerPortal, 2u);
    if (innerPortal != nullptr) {
      v->updatePortalPtr(innerPortal, 3u);
    }
  }

  // The new container bounds
  auto containerBounds =
      std::make_unique<CylinderVolumeBounds>(containerValues);

  // A z-binned volume link
  VariableVolumeLink vzLink(detail::VariableAxis(zBoundaries), binZ,
                            containerTransform.inverse());
  VariableVolumeLink ozLink(
      detail::VariableAxis(zBoundaries), binZ,
      outerPortal->surfaceRepresentation().transform(gctx).inverse());
  MultiplePortalLink oPortalMultiLink{oPortalLinks, std::move(ozLink)};

  // Update the outer portal (always exists)
  outerPortal->updatePortalLink(oPortalMultiLink, backward);
  // Update the inner portal (optional)
  if (innerPortal != nullptr) {
    VariableVolumeLink izLink(
        detail::VariableAxis(zBoundaries), binZ,
        innerPortal->surfaceRepresentation().transform(gctx).inverse());
    MultiplePortalLink iPortalMultiLink{iPortalLinks, std::move(izLink)};
    innerPortal->updatePortalLink(iPortalMultiLink, forward);
  }
  return DetectorVolume::makeShared(
      containerTransform, std::move(containerBounds),
      std::move(containerVolumes), std::move(vzLink), true, binZ, name);
}

std::shared_ptr<Acts::DetectorVolume>
Acts::CylindricalContainerBuilderPhi::build(const GeometryContext& gctx,
                                            const VolumeBuilder& volumeBuilder,
                                            const std::string& name) {
  return {};
}
