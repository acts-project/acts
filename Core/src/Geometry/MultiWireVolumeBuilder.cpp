// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/MultiWireVolumeBuilder.hpp"

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/IndexGrid.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

namespace Acts::Experimental {

MultiWireVolumeBuilder::MultiWireVolumeBuilder(
    const Config& config, std::unique_ptr<const Logger> logger)
    : m_config(config), m_logger(std::move(logger)) {
  if (m_config.mlSurfaces.empty()) {
    throw std::invalid_argument(
        "MultiWireStructureBuilder: No surfaces are given");
  }
}

std::unique_ptr<TrackingVolume> MultiWireVolumeBuilder::buildVolume() const {
  // Create the tracking volume

  ACTS_VERBOSE("Building a tracking volume with name "
               << m_config.name << " ,translation"
               << toString(m_config.transform.translation())
               << " and number of surfaces " << m_config.mlSurfaces.size());

  auto boundsType = m_config.bounds ? m_config.bounds->type()
                                    : VolumeBounds::BoundsType::eOther;
  if (!(boundsType == VolumeBounds::BoundsType::eTrapezoid ||
        boundsType == VolumeBounds::BoundsType::eCuboid ||
        boundsType == VolumeBounds::BoundsType::eDiamond)) {
    throw std::invalid_argument(
        "MultiWireStructureBuilder: Only trapezoid cuboid or diamond bounds "
        "are "
        "supported");
  }

  std::unique_ptr<TrackingVolume> trackingVolume =
      std::make_unique<TrackingVolume>(m_config.transform, m_config.bounds,
                                       m_config.name);

  // Add the surfaces to the tracking volume
  for (auto& surface : m_config.mlSurfaces) {
    trackingVolume->addSurface(surface);
  }

  return trackingVolume;
}

std::unique_ptr<Acts::NavigationPolicyFactory>
MultiWireVolumeBuilder::createNavigationPolicyFactory() const {
  if (m_config.binning.size() != 2u) {
    throw ::std::invalid_argument(
        "MultiWireStructureBuilder: Invalid binning provided");
  }
  auto [protoAxisA, expansionA] = m_config.binning.at(0);
  auto [protoAxisB, expansionB] = m_config.binning.at(1);

  // Create the grid from the axis
  const auto& iaxisA = protoAxisA.getAxis();
  const auto& iaxisB = protoAxisB.getAxis();
  // Binning needs to be equidistant
  if (iaxisA.getType() != AxisType::Equidistant ||
      iaxisB.getType() != AxisType::Equidistant) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: Binning axes need to be equidistant");
  }

  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisA(
      iaxisA.getBinEdges().front(), iaxisA.getBinEdges().back(),
      iaxisA.getNBins());

  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisB(
      iaxisB.getBinEdges().front(), iaxisB.getBinEdges().back(),
      iaxisB.getNBins());

  Grid<std::vector<std::size_t>, decltype(axisA), decltype(axisB)> grid(axisA,
                                                                        axisB);

  // The indexed grid to be filled from the navigation policy
  IndexGrid<decltype(grid)> indexedGrid(
      std::move(grid),
      {protoAxisA.getAxisDirection(), protoAxisB.getAxisDirection()},
      m_config.transform.inverse());

  TryAllNavigationPolicy::Config tryAllConfig;
  tryAllConfig.portals = true;
  tryAllConfig.sensitives = false;

  MultiLayerNavigationPolicy::Config navConfig;
  navConfig.binExpansion = {expansionA, expansionB};

  // Create the navigation policy factory
  std::unique_ptr<NavigationPolicyFactory> factory =
      NavigationPolicyFactory{}
          .add<TryAllNavigationPolicy>(tryAllConfig)
          .add<MultiLayerNavigationPolicy>(navConfig, indexedGrid)
          .asUniquePtr();

  return factory;
}

}  // namespace Acts::Experimental
