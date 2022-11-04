// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/NavigationDelegates.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceFactory.hpp"

#include <tuple>
#include <vector>

namespace {

/// @brief  Helper method to clip the surface to the given volume bounds
///
/// @param gctx the geometry context
/// @param volume the volume in question
/// @param surface the surface in question
///
/// @return an update transform and a new surface values container
std::tuple<Acts::Transform3, Acts::Surface::SurfaceType,
           Acts::SurfaceBounds::BoundsType, std::vector<Acts::ActsScalar>>
clipSurface(const Acts::GeometryContext& gctx,
            const Acts::Experimental::DetectorVolume& volume,
            const Acts::Surface& surface) {
  // primitve values
  const auto& vBounds = volume.volumeBounds();
  const auto vBoundValues = vBounds.values();
  std::vector<Acts::ActsScalar> bValues = surface.bounds().values();

  Acts::Surface::SurfaceType sType = surface.type();
  Acts::SurfaceBounds::BoundsType bType = surface.bounds().type();

  // Clip the cylinder
  if (vBounds.type() == Acts::VolumeBounds::BoundsType::eCylinder) {
    if (sType == Acts::Surface::SurfaceType::Cylinder) {
      bValues[Acts::CylinderBounds::BoundValues::eHalfLengthZ] =
          vBoundValues[Acts::CylinderVolumeBounds::BoundValues::eHalfLengthZ];
      Acts::Transform3 vTransform = volume.transform(gctx);
      return std::tie(vTransform, sType, bType, bValues);
    }
    if (surface.type() == Acts::Surface::SurfaceType::Disc) {
      bValues[Acts::RadialBounds::BoundValues::eMinR] =
          vBoundValues[Acts::CylinderVolumeBounds::BoundValues::eMinR];
      bValues[Acts::RadialBounds::BoundValues::eMaxR] =
          vBoundValues[Acts::CylinderVolumeBounds::BoundValues::eMaxR];
      Acts::Transform3 sTransform = surface.transform(gctx);
      return std::tie(sTransform, sType, bType, bValues);
    }
  }
  Acts::Transform3 sTransform = surface.transform(gctx);
  return std::tie(sTransform, sType, bType, bValues);
}

}  // namespace

std::vector<std::shared_ptr<Acts::Experimental::Portal>>
Acts::Experimental::DetrayJsonHelper::clipAndPickPortals(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::DetectorVolume& volume,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("JSON I/O:", logLevel));

  const auto& vPortals = volume.portals();
  std::vector<std::shared_ptr<Acts::Experimental::Portal>> clipped = {};

  ACTS_VERBOSE("Clip&pick portals of volume " << volume.name());

  for (auto [ip, p] : enumerate(vPortals)) {
    ACTS_VERBOSE("Processing portal " << ip);
    const auto& surface = p->surface();
    auto surfaceType = surface.type();
    // Check which direction is the keeper
    Acts::Vector3 outPointing(0., 0., 0.);
    Acts::Vector3 sNormal(0., 0., 0.);
    const Acts::Vector3 vCenter = volume.center(gctx);
    if (surfaceType == Acts::Surface::SurfaceType::Cylinder) {
      Acts::Vector3 vrPosition =
          vCenter + volume.volumeBounds().binningOffset(Acts::binR);
      Acts::Vector3 srPosition = surface.binningPosition(gctx, Acts::binR);
      sNormal = surface.normal(gctx, srPosition);
      outPointing = (srPosition - vrPosition);
    } else {
      const Acts::Vector3 sCenter = surface.center(gctx);
      sNormal = surface.normal(gctx, sCenter);
      outPointing = (sCenter - vCenter).normalized();
    }
    // Decide the keeping direction (for detray always outside)
    unsigned int outwards = (sNormal.dot(outPointing) > 0.) ? 1u : 0u;
    NavigationDirection oDir = (outwards == 1u) ? NavigationDirection::Forward
                                                : NavigationDirection::Backward;

    const auto& vLinks = p->volumeLinks();
    const auto& outLink = vLinks[outwards];

    // Determine the keeper
    const Acts::Experimental::IDelegateImpl* navImpl =
        outLink.implementation.get();

    // Clip surface to volume size and create portal
    auto [sTransform, sType, bType, bValues] =
        clipSurface(gctx, volume, surface);
    auto cSurface =
        SurfaceFactory::createSurface(sTransform, sType, bType, bValues);
    auto cPortal = Portal::makeShared(cSurface);
        clipped.push_back(cPortal);

    if (navImpl != nullptr) {
      /// @brief  Helper method to attach a single link
      auto attachSingleLink =
          [&](const Acts::Experimental::DetectorVolume& singleVolume) -> void {
        // Create a single outwards link
        auto nSingleLinkImpl =
            std::make_shared<detail::SingleLinkImpl>(singleVolume);
        DetectorVolumeLink nSingleLink;
        nSingleLink.connect<&detail::SingleLinkImpl::targetVolume>(
            nSingleLinkImpl.get());
        ManagedDetectorVolumeLink nManagedLink{std::move(nSingleLink),
                                               std::move(nSingleLinkImpl)};
        cPortal->updateVolumeLink(oDir, std::move(nManagedLink), {});
      };

      // Assign the correct volume links
      auto singleLink = dynamic_cast<const detail::SingleLinkImpl*>(navImpl);
      if (singleLink != nullptr) {
        ACTS_VERBOSE(
            "- Portal will have a single link implementation (native).");
        attachSingleLink(*singleLink->dVolume);
        continue;
      }
      // It is a multi link
      auto multiLink1D = dynamic_cast<const detail::MultiLink1DImpl*>(navImpl);
      if (multiLink1D != nullptr) {
        ActsScalar clippedMin = 0.;
        ActsScalar clippedMax = 0.;
        if (multiLink1D->bValue == binR) {
          clippedMin = bValues[RadialBounds::BoundValues::eMinR];
          clippedMax = bValues[RadialBounds::BoundValues::eMaxR];
        } else if (multiLink1D->bValue == binZ) {
          ActsScalar vZ = sTransform.translation().z();
          clippedMin = vZ - bValues[CylinderBounds::BoundValues::eHalfLengthZ];
          clippedMax = vZ + bValues[CylinderBounds::BoundValues::eHalfLengthZ];
        }
        // Clipp boundaries to the new range
        auto [clippedBoundaries, clippedVolumes] =
            clip(multiLink1D->cBoundaries, multiLink1D->dVolumes,
                 {clippedMin, clippedMax}, logLevel);
        // Check if after clipping a single link is enough
        if (clippedVolumes.size() == 1u) {
          ACTS_VERBOSE(
              "- Portal will have a single link implementation (clipped).");
          attachSingleLink(*clippedVolumes[0u]);
        } else {
          ACTS_VERBOSE(
              "- Portal will remain with a multi link implementation.");
          ACTS_VERBOSE("  " << clippedVolumes.size() << " of "
                            << multiLink1D->dVolumes.size() << " remain");
          ACTS_VERBOSE("  in clip range  ["
                       << clippedMin << ", " << clippedMax << "] out of ["
                       << multiLink1D->cBoundaries[0u] << ","
                       << multiLink1D->cBoundaries.back() << "]");
          // Create a single outwards link
          auto nMultiLinkImpl = std::make_shared<detail::MultiLink1DImpl>(
              clippedVolumes, clippedBoundaries, multiLink1D->bValue);
          DetectorVolumeLink nMultiLink;
          nMultiLink.connect<&detail::MultiLink1DImpl::targetVolume>(
              nMultiLinkImpl.get());
          ManagedDetectorVolumeLink nManagedLink{std::move(nMultiLink),
                                                 std::move(nMultiLinkImpl)};
          cPortal->updateVolumeLink(oDir, std::move(nManagedLink), {});
        }
        continue;
      }

      auto tMultiLink1D =
          dynamic_cast<const detail::TransformedMultiLink1DImpl*>(navImpl);
      if (tMultiLink1D != nullptr) {
        ActsScalar clippedMin = 0.;
        ActsScalar clippedMax = 0.;
        Transform3 tTransform = tMultiLink1D->transform;
        Vector3 tCenter = tTransform.translation();
        if (multiLink1D->bValue == binR) {
          clippedMin = bValues[RadialBounds::BoundValues::eMinR];
          clippedMax = bValues[RadialBounds::BoundValues::eMaxR];
        } else if (multiLink1D->bValue == binZ) {
          ActsScalar vZ = sTransform.translation().z();
          clippedMin = vZ - bValues[CylinderBounds::BoundValues::eHalfLengthZ];
          clippedMax = vZ + bValues[CylinderBounds::BoundValues::eHalfLengthZ];
        }

        // Clipp boundaries to the new range
        auto [clippedBoundaries, clippedVolumes] =
            clip(tMultiLink1D->multiLink.cBoundaries,
                 tMultiLink1D->multiLink.dVolumes, {clippedMin, clippedMax},
                 logLevel);
        // Check if after clipping a single link is enough
        if (clippedVolumes.size() == 1u) {
          ACTS_VERBOSE(
              "Portal will have a single link implementation (clipped).");
          attachSingleLink(*clippedVolumes[0u]);
        } else {
          ACTS_VERBOSE(
              "Portal will remain with a transformed multi link "
              "implementation.");
          // Create a single outwards link
          // auto MultiLinkImpl = std::make_shared<detail::MultiLink1DImpl>(
          //    clippedVolumes, clippedBoundaries, multiLink1D->bValue);
          // DetectorVolumeLink nMultiLink;
          // nMultiLink.connect<&detail::MultiLink1DImpl::targetVolume>(
          //    nMultiLinkImpl.get());
          // ManagedDetectorVolumeLink nManagedLink{std::move(nMultiLink),
          //                                       std::move(nMultiLinkImpl)};
          // cPortal->updateVolumeLink(oDir, std::move(nManagedLink), {});
        }
        continue;
      }
      // It is a multi link with transform or another one
      ACTS_WARNING("Not implemented yet, boundary will not be clipped.");
    }
  }

  return clipped;
}

std::tuple<std::vector<Acts::ActsScalar>,
           std::vector<const Acts::Experimental::DetectorVolume*>>
Acts::Experimental::DetrayJsonHelper::clip(
    const std::vector<Acts::ActsScalar>& vBoundaries,
    const std::vector<const Acts::Experimental::DetectorVolume*>& volumes,
    const std::array<Acts::ActsScalar, 2u>& clipRange,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("JSON I/O:", logLevel));

  std::vector<ActsScalar> clippedBoundaries = vBoundaries;
  std::vector<const Experimental::DetectorVolume*> clippedVolumes = volumes;

  if (clipRange[0u] > clippedBoundaries.front() or
      clipRange[1u] < clippedBoundaries.back()) {
    ACTS_DEBUG("Clip range smaller than original boundary range.");

    size_t lerase = 0;
    size_t herase = 0;
    // Erase borders and count
    for (auto it = clippedBoundaries.begin(); it != clippedBoundaries.end();) {
      if (*it <= clipRange[0u]) {
        ACTS_VERBOSE("- Removing " << *it
                                   << " as it is smaller (equal) than (to) "
                                   << clipRange[0u]);
        it = clippedBoundaries.erase(it);
        ++lerase;
      } else if (*it >= clipRange[1u]) {
        it = clippedBoundaries.erase(it);
        ACTS_VERBOSE("- Removing " << *it
                                   << " as it is bigger (equal) than (to) "
                                   << clipRange[1u]);
        ++herase;
      } else {
        ++it;
      }
    }
    // now erase the volumes not covered by the clip
    if (lerase > 1u) {
      ACTS_VERBOSE("- Erasing " << lerase - 1u
                                << " volumes at start of range.");
      clippedVolumes.erase(clippedVolumes.begin(),
                           clippedVolumes.begin() + lerase - 1u);
    }
    if (herase > 1u) {
      ACTS_VERBOSE("- Erasing " << herase - 1u << " volumes at end of range.");
      for (size_t ip = 1u; ip < herase; ++ip) {
        clippedVolumes.pop_back();
      }
    }
    // Set the new boundaries
    if (lerase > 0u) {
      clippedBoundaries.insert(clippedBoundaries.begin(), clipRange[0u]);
    }
    if (herase > 0u) {
      clippedBoundaries.push_back(clipRange[1u]);
    }

  } else {
    ACTS_DEBUG("Clip range bigger or equal to original boundary range.");
  }

  return std::tie(clippedBoundaries, clippedVolumes);
}
