// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelDetectorSurfaceFactory.hpp"

#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Utilities/Enumerate.hpp"

Acts::GeoModelDetectorSurfaceFactory::GeoModelDetectorSurfaceFactory(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

void Acts::GeoModelDetectorSurfaceFactory::construct(
    Cache& cache, const GeometryContext&, const GeoModelTree& geoModelTree,
    const Options& options) {
  if (geoModelTree.geoReader == nullptr) {
    throw std::invalid_argument("GeoModelTree has no GeoModelReader");
  }

  // The map cache
  std::vector<std::pair<GeometryIdentifier, std::shared_ptr<Surface>>> mapCache;

  // Loop over the queries and retrieve the published nodes
  for (const auto& q : options.queries) {
    ACTS_VERBOSE("Constructing detector elements for query " << q);
    auto qFPV =
        geoModelTree.geoReader->getPublishedNodes<std::string, GeoFullPhysVol*>(
            q);

    for (auto& [name, fpv] : qFPV) {
      // Convert using the list of converters
      for (const auto& converter : m_cfg.shapeConverters) {
        auto converted = converter->toSensitiveSurface(*fpv);
        if (converted.ok()) {
          // Add the element and surface to the cache
          cache.sensitiveSurfaces.push_back(converted.value());
        }
      }
    }
    ACTS_VERBOSE("Found " << qFPV.size()
                          << " full physical volumes matching the query.");
  }
  ACTS_DEBUG("Constructed "
             << cache.sensitiveSurfaces.size() << " sensitive elements and "
             << cache.passiveSurfaces.size() << " passive elements");
}
