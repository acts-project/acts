// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/MultiWireVolumeBuilder.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrapezoidPortalShell.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

// Constructor
Acts::MultiWireVolumeBuilder::MultiWireVolumeBuilder(
    const Config& config, std::unique_ptr<const Acts::Logger> logger)
    : m_config(config), m_logger(std::move(logger)) {
  if (m_config.mlSurfaces.empty()) {
    throw std::invalid_argument(
        "MultiWireStructureBuilder: No surfaces are given");
  } else if (m_config.binning.size() != 2u) {
    throw ::std::invalid_argument(
        "MultiWireStructureBuilder: Invalid binning provided");
  }
}

std::unique_ptr<Acts::TrackingVolume> Acts::MultiWireVolumeBuilder::buildVolume(
    Acts::GeometryContext& gctx) const {
  // Create the tracking volume

  ACTS_VERBOSE("Building a tracking volume with name "
               << m_config.name << " ,translation"
               << toString(m_config.transform.translation())
               << " and number of surfaces " << m_config.mlSurfaces.size());

  const auto& bounds =
      dynamic_pointer_cast<TrapezoidVolumeBounds>(m_config.bounds);
  if (bounds == nullptr) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: Invalid bounds - trapezoidal needed");
  }

  std::unique_ptr<Acts::TrackingVolume> trackingVolume =
      std::make_unique<Acts::TrackingVolume>(m_config.transform, bounds,
                                             m_config.name);

  SingleTrapezoidPortalShell portalShell(*trackingVolume);
  portalShell.applyToVolume();

  // Add the surfaces to the tracking volume
  for (auto& surface : m_config.mlSurfaces) {
    trackingVolume->addSurface(surface);
  }

  auto [protoAxisA, expansionA] = m_config.binning.at(0);
  auto [protoAxisB, expansionB] = m_config.binning.at(1);

  // Create the grid from the axis
  const auto& iaxisA = protoAxisA.getAxis();
  const auto& iaxisB = protoAxisB.getAxis();
  // Binning needs to be equidistant
  if (iaxisA.getType() != Acts::AxisType::Equidistant ||
      iaxisB.getType() != Acts::AxisType::Equidistant) {
    throw std::runtime_error(
        "MultiWireVolumeBuilder: Binning axes need to be equidistant");
  }

  Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound> axisA(
      iaxisA.getBinEdges().front(), iaxisA.getBinEdges().back(),
      iaxisA.getNBins());

  Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Bound> axisB(
      iaxisB.getBinEdges().front(), iaxisB.getBinEdges().back(),
      iaxisB.getNBins());

  Acts::Grid<std::vector<std::size_t>, decltype(axisA), decltype(axisB)> grid(
      axisA, axisB);

  // The indexed grid to be filled from the navigation policy
  Acts::Experimental::IndexedSurfacesNavigation<decltype(grid)> indexedGrid(
      std::move(grid),
      {protoAxisA.getAxisDirection(), protoAxisB.getAxisDirection()});

  // Use TryAll Navigation Policy for the portals and acceleration structure
  // with indexed surfaces for the sensitives

  Acts::TryAllNavigationPolicy::Config tryAllConfig;
  tryAllConfig.portals = true;
  tryAllConfig.sensitives = false;

  Acts::TryAllNavigationPolicy tryAllPolicy(gctx, *trackingVolume, *m_logger,
                                            tryAllConfig);

  // Configure the navigation policy with the binning for the grid for the
  // sensitive surfaces

  MultiLayerNavigationPolicy<decltype(indexedGrid)>::Config navConfig;
  navConfig.binExpansion = {expansionA, expansionB};
  navConfig.axis = {protoAxisA, protoAxisB};

  std::unique_ptr<NavigationPolicyFactory> factory =
      NavigationPolicyFactory::make()
          .add<TryAllNavigationPolicy>(tryAllConfig)
          .add<MultiLayerNavigationPolicy<decltype(indexedGrid)>>(navConfig,
                                                                  indexedGrid)
          .asUniquePtr();

  auto policyBase = factory->build(gctx, *trackingVolume, *m_logger);

  trackingVolume->setNavigationPolicy(std::move(policyBase));

  return trackingVolume;
}
