// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/GbtsLayerDescriptions.hpp"

#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "ActsPlugins/Json/GbtsConfigJsonConverter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace ActsExamples {

namespace {

/// Add a surface of the tracking geometry to the description of its GBTS
/// layer: the reference coordinates are summed and the bounds extended.
/// @param surface The surface to add
/// @param gctx The geometry context
/// @param layerMap The GBTS layer of the surfaces
/// @param inputVector The layer descriptions, one per GBTS layer seen so far
/// @param countVector The number of surfaces added to each layer description
/// @param fillModuleCsv Whether to write the module to ACTS_modules.csv
/// @param logger The logger a surface on none of the layers goes to
void addSurfaceToGbtsLayers(
    const Acts::Surface &surface, const Acts::GeometryContext &gctx,
    const GbtsLayerMap &layerMap,
    std::vector<Acts::Experimental::GbtsLayerDescription> &inputVector,
    std::vector<std::size_t> &countVector, bool fillModuleCsv,
    const Acts::Logger &logger) {
  Acts::GeometryIdentifier geoId = surface.geometryId();
  auto actsVolId = geoId.volume();
  auto actsLayId = geoId.layer();
  auto mod_id = geoId.sensitive();
  auto center = surface.center(gctx);

  // A polygonal surface gives its corners whatever the segment count
  // is, curved bounds an approximation of that resolution.
  const std::vector<Acts::Vector3> corners =
      surface.polyhedronRepresentation(gctx, 4u).vertices;

  float rc = 0.0;
  float minBound = std::numeric_limits<float>::infinity();
  float maxBound = -std::numeric_limits<float>::infinity();

  // convert to Gbts ID
  const auto find = layerMap.find(geoId);

  // a surface off the GBTS layers takes no part in the seeding
  if (find == layerMap.end()) {
    ACTS_DEBUG("No GBTS layer for volume: "
               << geoId.volume() << " Layer: " << geoId.layer()
               << " Surface: " << geoId.sensitive());
    return;  // skip this surface
  }

  const Acts::Experimental::GbtsExperimentLayerId gbtsId = find->layerId;

  // a variable that says if barrrel, 0 = barrel
  Acts::Experimental::GbtsLayerType barrelEc = find->type;

  if (barrelEc == Acts::Experimental::GbtsLayerType::Barrel) {
    rc = Acts::fastHypot(center.x(), center.y());  // barrel center in r
    // bounds of z
    for (const Acts::Vector3 &corner : corners) {
      minBound = std::min(minBound, static_cast<float>(corner.z()));
      maxBound = std::max(maxBound, static_cast<float>(corner.z()));
    }
  } else if (barrelEc == Acts::Experimental::GbtsLayerType::Endcap) {
    rc = center.z();  // not barrel center in Z
    // bounds of r
    for (const Acts::Vector3 &corner : corners) {
      const auto r =
          static_cast<float>(Acts::fastHypot(corner.x(), corner.y()));
      minBound = std::min(minBound, r);
      maxBound = std::max(maxBound, r);
    }
  } else {
    throw std::runtime_error("Invalid barrel/endcap assignment for GbtsLayer");
  }

  // tilt of a plane module in r-z, along which n_r dr + n_z dz = 0
  float slope = 0.f;
  if (surface.type() == Acts::Surface::Plane) {
    const Acts::Vector3 normal =
        surface.normal(gctx, center, Acts::Vector3::UnitZ());
    const double nr = normal.head<2>().dot(center.head<2>().normalized());
    const double nz = normal.z();
    slope = static_cast<float>(
        barrelEc == Acts::Experimental::GbtsLayerType::Barrel ? -nz / nr
                                                              : -nr / nz);
  }

  const auto currentIndex =
      find_if(inputVector.begin(), inputVector.end(),
              [gbtsId](auto n) { return n.id == gbtsId; });
  if (currentIndex != inputVector.end()) {  // not end so does exist
    const auto index = static_cast<std::size_t>(
        std::distance(inputVector.begin(), currentIndex));
    inputVector[index].refCoord += rc;
    inputVector[index].minBound =
        std::min(inputVector[index].minBound, minBound);
    inputVector[index].maxBound =
        std::max(inputVector[index].maxBound, maxBound);
    inputVector[index].slope += slope;
    countVector[index] += 1;  // increase count at the index

  } else {  // end so doesn't exists
    // make new if one with Gbts ID doesn't exist:
    inputVector.push_back(
        Acts::Experimental::GbtsLayerDescription{.id = gbtsId,
                                                 .type = barrelEc,
                                                 .technology = find->technology,
                                                 .refCoord = rc,
                                                 .minBound = minBound,
                                                 .maxBound = maxBound,
                                                 .slope = slope});
    // so the element exists and not divinding by 0
    countVector.push_back(1);
  }

  // add to file each time,
  // print to csv for each module, no repeats so dont need to make
  // map for averaging
  if (fillModuleCsv) {
    std::fstream fout;
    fout.open("ACTS_modules.csv", std::ios::out | std::ios::app);
    fout << actsVolId << ", "                        // vol
         << actsLayId << ", "                        // lay
         << mod_id << ", "                           // module
         << gbtsId << ","                            // Gbts id
         << center.z() << ", "                       // z
         << Acts::fastHypot(center.x(), center.y())  // r
         << "\n";
  }
}

}  // namespace

GbtsLayerMap makeGbtsLayerMap(const std::filesystem::path &layerMappingFile) {
  std::vector<GbtsLayerMap::InputElement> actsToGbtsMap;

  // one entry per surface of a layer, sensitive 0 for a whole geometry layer
  for (const Acts::Experimental::GbtsLayerConfig &layer :
       Acts::Experimental::readGbtsLayers(layerMappingFile)) {
    for (const Acts::GeometryIdentifier &surface : layer.surfaces) {
      const GbtsLayerInfo layerInfo{.layerId = layer.id,
                                    .type = layer.type,
                                    .technology = layer.technology};
      actsToGbtsMap.emplace_back(surface, layerInfo);
    }
  }

  return GbtsLayerMap(std::move(actsToGbtsMap));
}

std::vector<Acts::Experimental::GbtsLayerDescription> makeGbtsLayerDescriptions(
    const Acts::TrackingGeometry &trackingGeometry,
    const GbtsLayerMap &layerMap, const Acts::GeometryContext &gctx,
    const bool useLayerTilt, const bool fillModuleCsv,
    const Acts::Logger &logger) {
  std::vector<Acts::Experimental::GbtsLayerDescription> inputVector;
  std::vector<std::size_t> countVector;

  trackingGeometry.visitSurfaces([&](const Acts::Surface *surface) {
    addSurfaceToGbtsLayers(*surface, gctx, layerMap, inputVector, countVector,
                           fillModuleCsv, logger);
  });

  for (std::size_t i = 0; i < inputVector.size(); i++) {
    inputVector[i].refCoord = inputVector[i].refCoord / countVector[i];
    // mean module tilt; below 1e-3 it is rounding in the module placement
    const float slope = inputVector[i].slope / countVector[i];
    inputVector[i].slope =
        useLayerTilt && std::abs(slope) >= 1e-3f ? slope : 0.f;
  }

  return inputVector;
}

}  // namespace ActsExamples
