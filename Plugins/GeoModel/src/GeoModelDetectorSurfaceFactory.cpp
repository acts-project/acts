// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelDetectorSurfaceFactory.hpp"

#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/GeoModelSurfaceConverter.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"

using namespace Acts::detail;

Acts::GeoModelDetectorSurfaceFactory::GeoModelDetectorSurfaceFactory(
    std::unique_ptr<const Logger> mlogger)
    : m_logger(std::move(mlogger)) {}

void Acts::GeoModelDetectorSurfaceFactory::construct(
    Cache& cache, const GeometryContext& gctx, const GeoModelTree& geoModelTree,
    const Options& options) {
  if (geoModelTree.geoReader == nullptr) {
    throw std::invalid_argument("GeoModelTree has no GeoModelReader");
  }

  for (const auto& q : options.queries) {
    ACTS_VERBOSE("Constructing detector elements for query " << q);
    auto qFPV =
        geoModelTree.geoReader->getPublishedNodes<std::string, GeoFullPhysVol*>(
            q);

    for (auto& [name, fpv] : qFPV) {
      // Convert the full physical volume to a sensitive surface
      auto sensitive =
          GeoModelSurfaceConverter::convertToSensitiveSurface(*fpv);
      if (std::get<0>(sensitive) != nullptr) {
        // Add the surface to the cache
        cache.sensitiveSurfaces.push_back(sensitive);
      }
    }
    ACTS_VERBOSE("Found " << qFPV.size()
                          << " full physical volumes matching the query.");
  }
  ACTS_DEBUG("Constructed "
             << cache.sensitiveSurfaces.size() << " sensitive elements and "
             << cache.passiveSurfaces.size() << " passive elements");
}
