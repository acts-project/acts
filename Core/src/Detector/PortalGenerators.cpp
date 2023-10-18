// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/PortalGenerators.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <iterator>
#include <stdexcept>
#include <utility>

std::vector<std::shared_ptr<Acts::Experimental::Portal>>
Acts::Experimental::generatePortals(
    const Transform3& dTransform, const VolumeBounds& dBounds,
    const std::shared_ptr<DetectorVolume>& dVolume) noexcept(false) {
  if (dVolume == nullptr) {
    throw std::runtime_error("PortalsGenerator: no detector volume provided.");
  }
  // Get the oriented boundary surfaces
  auto orientedSurfaces = dBounds.orientedSurfaces(dTransform);

  // Create and fill the portal return vector
  std::vector<std::shared_ptr<Portal>> portals;
  for (auto [i, oSurface] : enumerate(orientedSurfaces)) {
    // Create a portal from the surface
    auto portal = Portal::makeShared(oSurface.first);
    // Create a shared link instance & delegate
    auto singleLinkImpl =
        std::make_unique<const SingleDetectorVolumeImpl>(dVolume.get());
    DetectorVolumeUpdator singleLink;
    singleLink.connect<&SingleDetectorVolumeImpl::update>(
        std::move(singleLinkImpl));
    // Update the volume link and the store
    portal->assignDetectorVolumeUpdator(oSurface.second, std::move(singleLink),
                                        {dVolume});
    // Portal is prepared
    portals.push_back(std::move(portal));
  }

  // The portals are returned
  return portals;
}

Acts::Experimental::PortalGenerator
Acts::Experimental::defaultPortalGenerator() {
  PortalGenerator pGenerator;
  pGenerator.connect<&generatePortals>();
  return pGenerator;
}

std::vector<std::shared_ptr<Acts::Experimental::Portal>>
Acts::Experimental::generatePortalsUpdateInternals(
    const Transform3& dTransform, const VolumeBounds& dBounds,
    const std::shared_ptr<DetectorVolume>& dVolume) noexcept(false) {
  if (dVolume == nullptr) {
    throw std::runtime_error(
        "generatePortalsUpdateInternals: no detector volume provided.");
  }

  // Setting link to the mother volume to all sub volumes of this volume
  for (auto& vPtr : dVolume->volumePtrs()) {
    for (auto& pPtr : vPtr->portalPtrs()) {
      // Creating a link to the mother
      auto motherLinkImpl =
          std::make_unique<const SingleDetectorVolumeImpl>(dVolume.get());
      DetectorVolumeUpdator motherLink;
      motherLink.connect<&SingleDetectorVolumeImpl::update>(
          std::move(motherLinkImpl));
      pPtr->assignDetectorVolumeUpdator(std::move(motherLink), {dVolume});
    }
  }
  // Return from the standard generator
  return generatePortals(dTransform, dBounds, dVolume);
}

Acts::Experimental::PortalGenerator
Acts::Experimental::defaultPortalAndSubPortalGenerator() {
  PortalGenerator pGenerator;
  pGenerator.connect<&generatePortalsUpdateInternals>();
  return pGenerator;
}
