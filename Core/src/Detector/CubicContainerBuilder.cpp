// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/CubicContainerBuilder.hpp"

#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Detector/detail/CubicDetectorHelper.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"

#include <algorithm>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace Acts {
namespace Experimental {
class DetectorVolume;
}  // namespace Experimental
}  // namespace Acts

Acts::Experimental::CubicContainerBuilder::CubicContainerBuilder(
    const Acts::Experimental::CubicContainerBuilder::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : IDetectorComponentBuilder(), m_cfg(cfg), m_logger(std::move(logger)) {
  // Check if builders are present
  if (m_cfg.builders.empty()) {
    throw std::invalid_argument(
        "CubicContainerBuilder: no sub builders provided.");
  }
  // Check if binning value is correctly chosen
  std::vector<BinningValue> allowed = {binX, binY, binZ};
  if (std::find(allowed.begin(), allowed.end(), m_cfg.binning) ==
      allowed.end()) {
    throw std::invalid_argument(
        "CubicContainerBuilder: invalid binning value.");
  }
}

Acts::Experimental::DetectorComponent
Acts::Experimental::CubicContainerBuilder::construct(
    const GeometryContext& gctx) const {
  // Return container object
  DetectorComponent::PortalContainer rContainer;
  bool atNavigationLevel = true;

  // Create the indivudal components, collect for both outcomes
  std::vector<DetectorComponent> components;
  ACTS_DEBUG("Building container from " << m_cfg.builders.size()
                                        << " components.");
  // Check through the component volumes - if every builder only
  // built exactly one volume, you are at pure navigation level
  // Collect the volumes
  std::vector<std::shared_ptr<DetectorVolume>> volumes;
  std::vector<DetectorComponent::PortalContainer> containers;
  std::vector<std::shared_ptr<DetectorVolume>> rootVolumes;
  // Run through the builders
  std::for_each(
      m_cfg.builders.begin(), m_cfg.builders.end(), [&](const auto& builder) {
        auto [cVolumes, cContainer, cRoots] = builder->construct(gctx);
        atNavigationLevel = (atNavigationLevel and cVolumes.size() == 1u);
        // Collect individual components, volumes, containers, roots
        volumes.insert(volumes.end(), cVolumes.begin(), cVolumes.end());
        containers.push_back(cContainer);
        rootVolumes.insert(rootVolumes.end(), cRoots.volumes.begin(),
                           cRoots.volumes.end());
      });
  // Navigation level detected, connect volumes (cleaner and faster than
  // connect containers)
  if (atNavigationLevel) {
    ACTS_VERBOSE(
        "Component volumes are at navigation level: connecting volumes.");
    // Connect volumes
    rContainer = detail::CubicDetectorHelper::connect(
        gctx, volumes, m_cfg.binning, {}, logger().level());
  } else {
    ACTS_VERBOSE("Components contain sub containers: connect containers.");
    // Connect containers
    rContainer = detail::CubicDetectorHelper::connect(
        gctx, containers, m_cfg.binning, {}, logger().level());
  }
  // Return the container
  return Acts::Experimental::DetectorComponent{
      {}, rContainer, RootDetectorVolumes{rootVolumes, tryRootVolumes()}};
}
