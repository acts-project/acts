// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/KDTreeTrackingGeometryBuilder.hpp"

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <functional>

Acts::KDTreeTrackingGeometryBuilder::KDTreeTrackingGeometryBuilder(
    const Acts::KDTreeTrackingGeometryBuilder::Config& cfg,
    std::unique_ptr<const Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  m_cfg.protoDetector.harmonize();
}

std::unique_ptr<const Acts::TrackingGeometry>
Acts::KDTreeTrackingGeometryBuilder::trackingGeometry(
    const GeometryContext& gctx) const {
  using MeasuredSurface =
      std::pair<std::array<ActsScalar, 2u>, std::shared_ptr<Surface>>;
  // Prepare all the surfaces
  std::vector<MeasuredSurface> surfacesMeasured;
  surfacesMeasured.reserve(m_cfg.surfaces.size());
  for (auto& s : m_cfg.surfaces) {
    auto ext = s->polyhedronRepresentation(gctx, 1u).extent();
    surfacesMeasured.push_back(MeasuredSurface{
        std::array<ActsScalar, 2u>{ext.medium(binZ), ext.medium(binR)}, s});
  }

  // Create the KDTree
  SurfaceKDT surfaceKDT(std::move(surfacesMeasured));

  // Walk through the proto builder
  auto protoWorld = m_cfg.protoDetector.worldVolume;

  KDTreeTrackingGeometryBuilder::Cache cCache;

  auto worldTrackingVolume =
      translateVolume(cCache, gctx, surfaceKDT, protoWorld);

  ACTS_INFO("Retrieved " << cCache.surfaceCounter
                         << " surfaces from the KDTTree.");

  // return the geometry to the service
  return std::make_unique<const Acts::TrackingGeometry>(worldTrackingVolume);
}

std::shared_ptr<Acts::TrackingVolume>
Acts::KDTreeTrackingGeometryBuilder::translateVolume(
    Cache& cCache, const GeometryContext& gctx, const SurfaceKDT& kdt,
    const ProtoVolume& ptVolume, const std::string& indent) const {
  ACTS_DEBUG(indent << "Processing ProtoVolume: " << ptVolume.name);
  std::vector<std::shared_ptr<const TrackingVolume>> translatedVolumes = {};

  if (not ptVolume.constituentVolumes.empty() and not ptVolume.layerContainer) {
    ACTS_VERBOSE(indent << "> container volume with "
                        << ptVolume.constituentVolumes.size()
                        << " constituents.");
    for (auto& cVolume : ptVolume.constituentVolumes) {
      translatedVolumes.push_back(translateVolume(
          cCache, gctx, kdt, cVolume, indent + m_cfg.hierarchyIndent));
    }
    // Make a container volume
    return m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
        gctx, translatedVolumes);
  }

  std::vector<std::shared_ptr<const Layer>> layers = {};
  if (not ptVolume.constituentVolumes.empty()) {
    ACTS_VERBOSE(indent << "> layer volume with "
                        << ptVolume.constituentVolumes.size()
                        << " layer volumes.");
    for (auto& plVolume : ptVolume.constituentVolumes) {
      layers.push_back(translateLayer(cCache, gctx, kdt, plVolume,
                                      indent + m_cfg.hierarchyIndent));
    }
    ACTS_VERBOSE("> translated into " << layers.size() << " layers.");
    // Create a volume
    auto rangeR = ptVolume.extent.range(Acts::binR);
    auto rangeZ = ptVolume.extent.range(Acts::binZ);
    // Create a new tracking volume with those layers
    return m_cfg.trackingVolumeHelper->createTrackingVolume(
        gctx, layers, {}, nullptr, rangeR.min(), rangeR.max(), rangeZ.min(),
        rangeZ.max(), ptVolume.name);
  }
  return nullptr;
}

/// @return a new tracking volume
std::shared_ptr<const Acts::Layer>
Acts::KDTreeTrackingGeometryBuilder::translateLayer(
    Cache& cCache, const GeometryContext& gctx, const SurfaceKDT& kdt,
    const ProtoVolume& plVolume, const std::string& indent) const {
  ACTS_DEBUG(indent + "Processing ProtoVolume: " << plVolume.name);
  // Try to pull from the kd tree
  RangeXD<2u, ActsScalar> zrRange;
  zrRange[0u] = plVolume.extent.range(Acts::binZ);
  zrRange[1u] = plVolume.extent.range(Acts::binR);

  auto layerSurfaces = kdt.rangeSearchWithKey(zrRange);
  ACTS_VERBOSE(indent + ">> looking z/r range = " << zrRange.toString());
  ACTS_VERBOSE(indent + ">> found " << layerSurfaces.size()
                                    << " surfaces in the KDTree.");
  cCache.surfaceCounter += layerSurfaces.size();

  // Make a const collection out of the surfaces
  std::vector<std::shared_ptr<const Surface>> cLayerSurfaces;
  cLayerSurfaces.reserve(layerSurfaces.size());
  for (const auto& s : layerSurfaces) {
    cLayerSurfaces.push_back(s.second);
  }

  if (plVolume.layerType == Acts::Surface::SurfaceType::Cylinder) {
    ACTS_VERBOSE(indent + ">> creating cylinder layer.");
    return m_cfg.layerCreator->cylinderLayer(
        gctx, cLayerSurfaces, Acts::equidistant, Acts::equidistant);
  } else if (plVolume.layerType == Acts::Surface::SurfaceType::Disc) {
    ACTS_VERBOSE(indent + ">> creating disc layer.");
    return m_cfg.layerCreator->discLayer(gctx, cLayerSurfaces,
                                         Acts::equidistant, Acts::equidistant);
  }

  return nullptr;
}