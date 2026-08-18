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

namespace Acts {

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

  std::unique_ptr<TrackingVolume> trackingVolume{};

  if (m_config.alignablePlacement == nullptr) {
    trackingVolume = std::make_unique<TrackingVolume>(
        m_config.transform, m_config.bounds, m_config.name);
  } else {
    trackingVolume = std::make_unique<TrackingVolume>(
        *m_config.alignablePlacement, m_config.bounds, m_config.name);
  }

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
  auto [axisSpecA, expansionA] = m_config.binning.at(0);
  auto [axisSpecB, expansionB] = m_config.binning.at(1);

  // Binning needs to be fully specified and equidistant, with a direction
  if (axisSpecA.isDeferred() || axisSpecB.isDeferred()) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: Binning axes need a fully specified range");
  }
  if (!axisSpecA.isEquidistant() || !axisSpecB.isEquidistant()) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: Binning axes need to be equidistant");
  }
  if (!axisSpecA.direction().has_value() ||
      !axisSpecB.direction().has_value()) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: Binning axes need a direction");
  }

  // Create the grid from the axis
  const auto& paramsA = axisSpecA.asEquidistant();
  const auto& paramsB = axisSpecB.asEquidistant();

  // isDeferred() above guarantees the ranges are set
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisA(
      *paramsA.min, *paramsA.max, paramsA.nBins);

  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisB(
      *paramsB.min, *paramsB.max, paramsB.nBins);

  Grid<std::vector<std::size_t>, decltype(axisA), decltype(axisB)> grid(axisA,
                                                                        axisB);

  // The indexed grid to be filled from the navigation policy
  const auto* placement = m_config.alignablePlacement;
  auto indexedGrid =
      placement == nullptr
          ? IndexGrid<decltype(grid)>{std::move(grid),
                                      {*axisSpecA.direction(),
                                       *axisSpecB.direction()},
                                      m_config.transform.inverse()}
          : IndexGrid<decltype(grid)>{
                std::move(grid),
                {*axisSpecA.direction(), *axisSpecB.direction()},
                [placement](const GeometryContext& gctx) -> const Transform3& {
                  return placement->globalToLocalTransform(gctx);
                }};

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

}  // namespace Acts
