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
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"

#include <stdexcept>

std::vector<std::shared_ptr<Acts::Experimental::Portal>>
Acts::Experimental::DefaultPortalGenerator::generate(
    const Transform3& dTransform, const VolumeBounds& dBounds,
    const std::shared_ptr<DetectorVolume>& dVolume) const {
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
    auto singleLink =
        makeDetectorVolumeFinder<const SingleDetectorVolume>(dVolume.get());
    // Update the volume link and the store
    portal->assignDetectorVolumeFinder(oSurface.second, std::move(singleLink),
                                       {dVolume});
    // Portal is prepared
    portals.push_back(std::move(portal));
  }

  // The portals are returned
  return portals;
}

std::vector<std::shared_ptr<Acts::Experimental::Portal>>
Acts::Experimental::PortalAndSubPortalGenerator::generate(
    const Transform3& dTransform, const VolumeBounds& dBounds,
    const std::shared_ptr<DetectorVolume>& dVolume) const {
  if (dVolume == nullptr) {
    throw std::runtime_error(
        "generatePortalsUpdateInternals: no detector volume provided.");
  }

  // Setting link to the mother volume to all sub volumes of this volume
  for (auto& vPtr : dVolume->volumePtrs()) {
    for (auto& pPtr : vPtr->portalPtrs()) {
      // Creating a link to the mother
      auto motherLink =
          makeDetectorVolumeFinder<const SingleDetectorVolume>(dVolume.get());
      pPtr->assignDetectorVolumeFinder(std::move(motherLink), {dVolume});
    }
  }

  // Return from the standard generator
  return defaultGenerator.generate(dTransform, dBounds, dVolume);
}
