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

#include <GeoModelKernel/GeoShapeUnion.h>
#include <GeoModelKernel/GeoShapeShift.h>

namespace {
  std::string gname(const GeoShapeShift &gshift);
  std::string gname(const GeoShapeUnion &gunion);
  std::string gname(const GeoShape &gshape);

  std::string gname(const GeoShapeShift &gshift){
    return "Shift[" + gname(*gshift.getOp()) + "]";
  }
  std::string gname(const GeoShapeUnion &gunion){
    return "Union[" + gname(*gunion.getOpA()) + ", " + gname(*gunion.getOpB()) + "]";
  }
  std::string gname(const GeoShape &gshape){
    if( auto ptr = dynamic_cast<const GeoShapeUnion *>(&gshape); ptr != nullptr){
      return gname(*ptr);
    }
    if( auto ptr = dynamic_cast<const GeoShapeShift *>(&gshape); ptr != nullptr) {
      return gname(*ptr);
    }
    return gshape.type();
  }
}

Acts::GeoModelDetectorSurfaceFactory::GeoModelDetectorSurfaceFactory(
    const Config& cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

void Acts::GeoModelDetectorSurfaceFactory::construct(
    Cache& cache, const GeometryContext&, const GeoModelTree& geoModelTree,
    const Options& options) {
  if (geoModelTree.geoReader == nullptr) {
    throw std::invalid_argument("GeoModelTree has no GeoModelReader");
  }

  for (const auto& q : options.queries) {
    ACTS_VERBOSE("Constructing detector elements for query " << q);
    auto qFPV =
        geoModelTree.geoReader->getPublishedNodes<std::string, GeoFullPhysVol*>(
            q);

    auto matches = [&](const std::string &str) {
      if( m_cfg.nameList.empty() ) {
        return true;
      }
      return std::any_of(m_cfg.nameList.begin(), m_cfg.nameList.end(), [&](const auto &n){ return str.find(n) != std::string::npos; });
    };

    for (auto& [name, fpv] : qFPV) {
      const auto &vname = fpv->getLogVol()->getName();
      const auto &shape = *fpv->getLogVol()->getShape();

      // Mask
      if( !matches(name) ) {
        continue;
      }
      // Convert using the list of converters
      bool success = false;
      for (const auto& converter : m_cfg.shapeConverters) {
        auto converted = converter->toSensitiveSurface(*fpv);
        if (converted.ok()) {
          // Add the element and surface to the cache
          cache.sensitiveSurfaces.push_back(converted.value());
          success = true;
          ACTS_VERBOSE("successfully converted " << name << " (" << vname << " / " << gname(shape) << ")");
          break;
        }
      }

      if (!success) {
        ACTS_DEBUG(name << " (" << vname << " / " << gname(shape) << ") could not be converted by any converter");
      }
    }
    ACTS_VERBOSE("Found " << qFPV.size()
                          << " full physical volumes matching the query.");
  }
  ACTS_DEBUG("Constructed "
             << cache.sensitiveSurfaces.size() << " sensitive elements and "
             << cache.passiveSurfaces.size() << " passive elements");
}
