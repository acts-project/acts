// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/detail/PortalGenerators.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <tuple>
#include <vector>

namespace Acts {
namespace Experimental {
namespace detail {

/// @brief Calls the portal generation and adds registration to sub portals
///
/// This code is split off the PortalGenerator code in order to allow
/// unit testing of the portal generation wihtout detector volume construction
///
/// @param dTransform the context-resolved transform of the detector volume
/// @param dBounds the detecor volume bounds
/// @param dVolume the reference to the detector volume which generates this volume
///
/// @return a vector of newly created portals with registered inside volume
inline static std::vector<std::shared_ptr<Portal>> portalsAndSubPortals(
    const Transform3& dTransform, const VolumeBounds& dBounds,
    std::shared_ptr<DetectorVolume> dVolume) {
  // Setting link to the mother volume to all sub volumes of this volume
  for (auto vPtr : dVolume->volumePtrs()) {
    for (auto pPtr : vPtr->portalPtrs()) {
      // Creatint a link to the mother
      auto motherLinkImpl =
          std::make_shared<SingleDetectorVolumeImpl>(dVolume.get());
      DetectorVolumeUpdator motherLink;
      motherLink.connect<&SingleDetectorVolumeImpl::update>(
          motherLinkImpl.get());
      // Set it ot the portal
      ManagedDetectorVolumeUpdator managedMotherLink{std::move(motherLink),
                                                     std::move(motherLinkImpl)};
      pPtr->assignDetectorVolumeUpdator(std::move(managedMotherLink),
                                        {dVolume});
    }
  }
  // Return from the standard generator
  return portals(dTransform, dBounds, dVolume);
}

/// Create a default portal generator that connects to the
/// static method.
///
inline static Delegate<std::vector<std::shared_ptr<Portal>>(
    const Transform3&, const VolumeBounds&, std::shared_ptr<DetectorVolume>)>
defaultPortalAndSubPortalGenerator() {
  Delegate<std::vector<std::shared_ptr<Portal>>(
      const Transform3&, const VolumeBounds&, std::shared_ptr<DetectorVolume>)>
      pGenerator;
  pGenerator.connect<&portalsAndSubPortals>();
  return pGenerator;
}

}  // namespace detail
}  // namespace Experimental
}  // namespace Acts